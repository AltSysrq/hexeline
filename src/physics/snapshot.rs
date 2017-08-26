//-
// Copyright (c) 2017 Jason Lingle
//
// Permission to  use, copy,  modify, and/or distribute  this software  for any
// purpose  with or  without fee  is hereby  granted, provided  that the  above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE  IS PROVIDED "AS  IS" AND  THE AUTHOR DISCLAIMS  ALL WARRANTIES
// WITH  REGARD   TO  THIS  SOFTWARE   INCLUDING  ALL  IMPLIED   WARRANTIES  OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT  SHALL THE AUTHOR BE LIABLE FOR ANY
// SPECIAL,  DIRECT,   INDIRECT,  OR  CONSEQUENTIAL  DAMAGES   OR  ANY  DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
// OF  CONTRACT, NEGLIGENCE  OR OTHER  TORTIOUS ACTION,  ARISING OUT  OF OR  IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

/**

Logic for handling snapshots of the game state.

A "snapshot" is simply the state of the game at a particular frame. The concept
exists since we need to be able to revert to a previous snapshot if an event
comes in after the frame that event applies to.

Sequential snapshots will attempt to share heap memory, primarily to reduce
the memory bandwidth that would be consumed by needing to copy the whole heap
every frame.

In order to support parallel processing of the "update objects" and "collision
detection" stages of the update pipeline, the snapshot system is also tightly
coupled to some other components that would otherwise be separate.

*/

use std::mem;
use std::sync::Arc;

use simd::*;

use super::bounding_rhombus::BoundingRhombus;
use super::common_object::CommonObject;

/// Manages the raw memory allocation behind an arena.
///
/// This is essentially a `Vec` in disguise.
///
/// We use a dynamic memory expansion strategy by at least doubling `capacity`
/// when space is exhausted. This allows us to avoid allocating the full 1MB in
/// cases where it can be avoided. Note that these expansions do not entail an
/// evacuation, and end up bringing obsolete data long, since the arena has no
/// way to determine what data is live, much less update pointers to it.
///
/// Note that memory expansion does _not_ mutate the `ArenaMemory`; we instead
/// create an entirely new one, which allows pointers into the old one to
/// remain valid.
#[derive(Debug)]
struct ArenaMemory {
    base: *mut i32x4,
    capacity: usize,
}

unsafe impl Sync for ArenaMemory { }
unsafe impl Send for ArenaMemory { }

impl ArenaMemory {
    fn with_capacity(cap: usize) -> Self {
        assert!(cap <= 65536);

        let mut vec = Vec::with_capacity(cap);
        let ret = ArenaMemory {
            base: vec.as_mut_ptr(),
            capacity: vec.capacity(),
        };
        mem::forget(vec);
        ret
    }
}

impl Drop for ArenaMemory {
    fn drop(&mut self) {
        mem::drop(unsafe { Vec::from_raw_parts(self.base, 0, self.capacity) });
    }
}

pub struct Snapshot {
    /// The objects in this snapshot, in canonical order.
    ///
    /// Canonical order is essentially "the way the system happens to lay the
    /// objects out", but does need to be consistent across all members in a
    /// networked game. The canonical order of a snapshot is largely based on
    /// the state of the previous snapshot.
    ///
    /// - All objects in both the new and old snapshots are in the same
    /// relative order to each other.
    ///
    /// - Objects emitted by other objects are ordered immediately after their
    /// parent object, in the order they were emitted.
    ///
    /// - Objects spontaneously emitted are added to the end in the order they
    /// were emitted.
    ///
    /// An exception is that on particular frames, the array is re-sorted
    /// according to modified Hilbert order.
    objects: Vec<CommonObject>,
    /// A tree of bounding rhombi used for collision detection.
    ///
    /// This is exactly parallel to `objects`.
    ///
    /// Each entry is the bounding rhombus of the corresponding entry in
    /// `objects` unioned with its two children in this tree, if any. The tree
    /// is stored as a _reversed_ implicit binary tree; that is, children
    /// precede their parents. The tree can be traversed similarly to a normal
    /// implicit binary tree, but operating on negative offsets from a logical
    /// root. The logical root of a node N is M-2, where M is the least power
    /// of two greater than or equal to N+2. The direct children of the root
    /// are found at indices M-3 and M/2; in general, the children of a node D
    /// at depth K are found at D-1 and D-M/2**K, starting with K=1 at the
    /// root. A node is a leaf if M/2**K is zero.
    ///
    /// Note that since the tree is inverted, it is not necessarily complete.
    /// I.e., a node can have a sibling but no parent, or even multiple other
    /// nodes that should be descendents of its sibling even though that
    /// sibling does not exist.
    rhombus_tree: Vec<BoundingRhombus>,
    /// The arena held by this `Snapshot`.
    arena: Arc<ArenaMemory>,
    /// The number of entries in `arena` which are reachable from this
    /// `Snapshot`.
    ///
    /// Entries beyond this index are volatile and accessing them is undefined
    /// behaviour.
    arena_size: usize,
}
