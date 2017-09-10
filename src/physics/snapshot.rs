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

use std::cmp::min;
use std::collections::{BTreeSet, HashMap};
use std::intrinsics::prefetch_read_data;
use std::mem;
use std::ptr;
use std::slice;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::atomic::hint_core_should_pause;

use arrayvec::ArrayVec;
use crossbeam::sync::ArcCell;
use odds::{get_unchecked, get_unchecked_mut, slice_unchecked};
use simd::*;
use simdext::*;

use super::bounding_rhombus::{BoundingRhombus, InvertedBoundingRhombus};
use super::common_object::CommonObject;
use super::extended_object::ExtendedObject;
use super::coords::*;
use super::composite::*;
use super::event::ObjectEvent;
use super::id_free_list::*;
use super::implicit_tree::*;
use super::hilbert_sort::xy_to_hilbert;
use super::units::*;
use super::xform::Affine2dH;

/// Type passed to `ExetndedObject::update()` to allow spawning new objects.
pub type SpawnList = ArrayVec<[CommonObject;32]>;
/// A set of object events which occur in a single frame.
pub type FrameObjectEvents = HashMap<u16, BTreeSet<ObjectEvent>>;

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

    #[inline]
    unsafe fn decompress_ptr(&self, ptr: u16, len: u8) -> &[i32x4] {
        let ptr = ptr as usize;
        let len = len as usize;
        debug_assert!(ptr + len <= self.capacity);
        slice::from_raw_parts((self.base as *const i32x4)
                              .offset(ptr as isize), len)
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
    /// There are two padding entries at the beginning to allow eliding bounds
    /// checks in particular cases.
    ///
    /// This is exactly parallel to `objects`.
    ///
    /// Each entry is the bounding rhombus of the corresponding entry in
    /// `objects` unioned with its two children in this tree, if any. The tree
    /// is stored as a _reversed_ implicit binary tree, as per
    /// `ImplicitTreeNavigator`.
    ///
    /// Note that since the tree is inverted, it is not necessarily complete.
    /// I.e., a node can have a sibling but no parent, or even multiple other
    /// nodes that should be descendents of its sibling even though that
    /// sibling does not exist.
    rhombus_tree: Vec<BoundingRhombus>,
    /// An `ImplicitTreeNavigator` at the top root of `rhombus_tree`.
    rhombus_tree_navigator: ImplicitTreeNavigator,
    /// The arena held by this `Snapshot`.
    arena: Arc<ArenaMemory>,
    /// The number of entries in `arena` which are reachable from this
    /// `Snapshot`.
    ///
    /// Entries beyond this index are volatile and accessing them is undefined
    /// behaviour.
    arena_size: usize,
    /// Collisions detected on the previous frame, sorted. Each collision adds
    /// two entries, one for each object.
    collisions: Vec<ImpendingCollision>,

    free_ids: IdFreeList,

    /// The frame number of this snapshot.
    tick: u32,
}

/// A non-commutative collision detected between two objects.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ImpendingCollision {
    /// The index of the object affected by the collision.
    subject: u16,
    /// The index of the other object involved in the collision.
    other: u16,
}

impl Snapshot {
    /// Decompress a pointer into this snapshot's arena.
    ///
    /// ## Unsafety
    ///
    /// Bounds checking is only performed in debug builds.
    #[inline(always)]
    pub unsafe fn decompress_ptr(&self, ptr: u16, len: u8) -> &[i32x4] {
        debug_assert!(ptr as usize + len as usize <= self.arena_size);
        self.arena.decompress_ptr(ptr, len)
    }

    /// Like `decompress_ptr`, but takes input as a `CommonObject`.
    #[inline(always)]
    pub unsafe fn decompress_obj(&self, obj: CommonObject) -> &[i32x4] {
        self.decompress_ptr(obj.extended_data(), obj.data_dst_size())
    }

    #[inline(always)]
    pub fn tick(&self) -> u32 { self.tick }
}

#[inline]
fn find_intersecting<'a>(navigator: ImplicitTreeNavigator,
                         objects: &'a [CommonObject],
                         rhombus_tree: &'a [BoundingRhombus],
                         bound: InvertedBoundingRhombus)
                         -> impl Iterator<Item = (u32, CommonObject)> + 'a {
    navigator.traverse_rp(
        move |ix| rhombus_tree.get(ix as usize + 2).map_or(
            true, |rhombus| rhombus.overlaps(bound)),
        move |ix| unsafe {
            let ix = ix as isize;
            prefetch_read_data(rhombus_tree.as_ptr().offset(ix+2), 3);
            prefetch_read_data(objects.as_ptr().offset(ix), 2);
        })
        .filter_map(move |ix| objects.get(ix as usize).cloned()
                    .map(move |o| (ix, o)))
        .filter(move |&(_, object)|
                unlikely!(object.bounding_rhombus().overlaps(bound)))
}

/// Centralised state for building a new `Snapshot`.
///
/// Building a snapshot takes place in a pipeline of two stages, referred to as
/// "update" and "collision", which may be run in parallel or sequentially.
pub struct SnapshotWriter {
    /// The previous snapshot used as a basis for this one.
    prev_snapshot: Arc<Snapshot>,

    /// Shared mutable state.
    ///
    /// The vectors in `objects` do not have any particular length. They are
    /// not really used as vectors here, but `Vec` is used to handle memory
    /// management. Manipulation of the contents is done with unsafe pointer
    /// accesses. The vectors are not resized in-place, but are reallocated
    /// entirely so that the memory held by the collision stage remains valid
    ///
    /// `num_objects` is a monotonically increasing value indicating how many
    /// values have been written to `objects` by the update stage. Reading
    /// `num_objects` and then `objects` allows the collision stage to get a
    /// slice of both values which are initialised and stable.
    ///
    /// Objects which take part in collisions need to have their wakeup_counter
    /// set to 255. It might seem like the collision stage could do this itself
    /// by writing to indices below `num_objects`, but this is not sound as the
    /// update stage needs to be able to clone the array when it reallocates.
    objects: ArcCell<(Vec<CommonObject>, Vec<BoundingRhombus>)>,
    num_objects: AtomicUsize,
    /// The current arena memory. Similarly to `objects` et al, this is
    /// actively appended to by the update stage while the collision stage may
    /// be reading from it. There is know way for the collision stage to know
    /// the actual arena size; it must simply trust that any references in the
    /// initialised portion of `objects` are in-bounds.
    arena: ArcCell<ArenaMemory>,
    /// The amount of `arena` already in use by other snapshots.
    init_arena_size: usize,

    /// Flags indicating whether the respective pipeline stages have been
    /// started. Used for sanity checks.
    update_started: AtomicBool,
    collision_started: AtomicBool,

    /// Flag indicating whether the update stage has completed.
    ///
    /// When all available objects have been processed, the collision stage
    /// spin waits until either more objects are available or the update stage
    /// has completed.
    update_done: AtomicBool,
}

/// The output from the update stage of the snapshot pipeline.
pub struct SnapshotUpdateResult {
    num_objects: usize,
    arena_size: usize,
    free_ids: IdFreeList,
}

/// The output from the collision stage of the snapshot pipeline.
pub struct SnapshotCollisionResult {
    collisions: Vec<ImpendingCollision>,
}

impl SnapshotWriter {
    /// Create a new `SnapshotWriter` at frame 1, with a `prev_snapshot` that
    /// is completely empty.
    pub fn empty(object_capacity: usize, arena_capacity: usize) -> Self {
        let padding = BoundingRhombus::around(Vhs(0, 0), 0);
        let arena = Arc::new(ArenaMemory::with_capacity(arena_capacity));
        let free_ids = IdFreeListBuilder::new().build();

        SnapshotWriter::wrap(Arc::new(Snapshot {
            objects: Vec::new(),
            rhombus_tree: vec![padding, padding],
            rhombus_tree_navigator: ImplicitTreeNavigator::root_of(0),
            arena: arena.clone(),
            arena_size: 0,
            collisions: Vec::new(),
            free_ids: free_ids.clone(),
            tick: 0,
        }), object_capacity)
    }

    /// Create a writer for the frame after `prev_snapshot`.
    ///
    /// This must only be used if no snapshots in the same chain exist with a
    /// later tick than `prev_snapshot`.
    fn wrap(prev_snapshot: Arc<Snapshot>, object_cap: usize) -> Self {
        SnapshotWriter {
            objects: ArcCell::new(Arc::new(
                (Vec::with_capacity(object_cap),
                 Vec::with_capacity(2 + object_cap)))),
            num_objects: AtomicUsize::new(0),
            arena: ArcCell::new(prev_snapshot.arena.clone()),
            init_arena_size: prev_snapshot.arena_size,
            update_started: AtomicBool::new(false),
            collision_started: AtomicBool::new(false),
            update_done: AtomicBool::new(false),

            prev_snapshot: prev_snapshot,
        }
    }

    /// Create an independent snapshot writer building the frame after
    /// `prev_snapshot`.
    ///
    /// It is permissible for there to be snapshots in the chain newer than
    /// `prev_snapshot`, but note that this call is somewhat expensive as a
    /// result.
    pub fn fork(prev_snapshot: Arc<Snapshot>) -> Self {
        let arena = ArenaMemory::with_capacity(
            prev_snapshot.arena.capacity);
        unsafe {
            ptr::copy_nonoverlapping(
                prev_snapshot.arena.base,
                arena.base,
                prev_snapshot.arena_size);
        }
        let num_objects = prev_snapshot.objects.len();

        SnapshotWriter {
            objects: ArcCell::new(Arc::new(
                (Vec::with_capacity(num_objects + 256),
                 Vec::with_capacity(2 + num_objects + 256)))),
            num_objects: AtomicUsize::new(0),
            arena: ArcCell::new(Arc::new(arena)),
            init_arena_size: prev_snapshot.arena_size,
            update_started: AtomicBool::new(false),
            collision_started: AtomicBool::new(false),
            update_done: AtomicBool::new(false),

            prev_snapshot: prev_snapshot,
        }
    }

    /// Returns the snapshot prior to the one being written.
    pub fn prev_snapshot(&self) -> &Arc<Snapshot> {
        &self.prev_snapshot
    }

    /// Fully run the update stage of the snapshot pipeline.
    ///
    /// `extra` is invoked after the body of the pipeline has completed. This
    /// callback can be used to spawn extra objects. It is passed a function
    /// that it can call to spawn such objects.
    ///
    /// Returns the partial result of this pipeline stage.
    #[inline(never)]
    pub fn run_update(&self, events: &FrameObjectEvents,
                      // Use dynamic dispatch since this function will be quite
                      // large, so we don't want to instantiate it multiple
                      // times. Also, performance of the callback is
                      // irrelevant.
                      extra: &mut for<'b> FnMut (
                          &mut SnapshotUpdatePipeline<'b>,
                          fn (&mut SnapshotUpdatePipeline<'b>, CommonObject)))
                      -> SnapshotUpdateResult {
        assert!(!self.update_started.swap(true, Ordering::Relaxed));

        let mut pipeline = SnapshotUpdatePipeline::new(self, events);
        pipeline.run();
        extra(&mut pipeline, SnapshotUpdatePipeline::add_object);

        // Done, make sure the collision stage is made aware
        self.num_objects.store(pipeline.result.num_objects, Ordering::Release);
        self.update_done.store(true, Ordering::Release);

        pipeline.result
    }
}

pub struct SnapshotUpdatePipeline<'a> {
    writer: &'a SnapshotWriter,
    events: &'a FrameObjectEvents,

    // The current working memory.
    objects: Arc<(Vec<CommonObject>, Vec<BoundingRhombus>)>,
    arena: Arc<ArenaMemory>,

    // The final number of objects we currently expect to have.
    reserved: usize,
    // The number of objects we've allocated space for.
    capacity: usize,

    tree_builder: ImplicitTreeBuilder,
    prev_bounding_rhombus: BoundingRhombus,

    collisions_in: &'a [ImpendingCollision],

    // Accumulator for the final result. `num_objects` is used to track the
    // running index of the current output object.
    result: SnapshotUpdateResult,
}

impl<'a> SnapshotUpdatePipeline<'a> {
    fn new(writer: &'a SnapshotWriter, events: &'a FrameObjectEvents) -> Self {
        let objects = writer.objects.get();
        let capacity = objects.0.capacity();
        assert_eq!(capacity + 2, objects.1.capacity());
        assert!(writer.prev_snapshot.objects.len() <= capacity);

        SnapshotUpdatePipeline {
            writer, events, objects, capacity,
            reserved: writer.prev_snapshot.objects.len(),
            arena: writer.arena.get(),
            tree_builder: ImplicitTreeBuilder::new(),
            prev_bounding_rhombus: BoundingRhombus::around(Vhs(0, 0), 0),
            collisions_in: &writer.prev_snapshot.collisions,
            result: SnapshotUpdateResult {
                num_objects: 0,
                arena_size: writer.init_arena_size,
                free_ids: writer.prev_snapshot.free_ids.clone(),
            },
        }
    }

    fn run(&mut self) {
        let tick_mod_4 = (self.writer.prev_snapshot.tick + 1) & 3;

        for (src_index, mut object) in
            self.writer.prev_snapshot.objects.iter().cloned().enumerate()
        {
            object = object.tick(tick_mod_4 as i32);
            if unlikely!(0 == object.wakeup_counter()) {
                self.update_extended(src_index, object);
            } else {
                self.write_object(object);
            }

            // Periodically update the object count to allow the collision
            // stage to continue.
            //
            // We only do this every 64 objects, and additionally offset the
            // threshold by 64 objects, not only to reduce the global bus
            // events on architectures with weaker memory ordering, but also to
            // ensure that the collision stage doesn't actively read the same
            // cache lines we're still writing.
            if 0 == src_index % 64 && self.result.num_objects > 64 {
                self.writer.num_objects.store(
                    self.result.num_objects - 64,
                    Ordering::Release);
            }
        }

        assert_eq!(self.reserved, self.result.num_objects);
    }

    fn update_extended(&mut self, src_ix: usize, object: CommonObject) {
        let collision_split_point = self.collisions_in.iter().enumerate()
            .take_while(|&(_, c)| c.subject == src_ix as u16)
            .map(|(ix, _)| ix)
            .last()
            .unwrap_or(0);
        let (this_collisions, other_collisions) =
            self.collisions_in.split_at(collision_split_point);
        self.collisions_in = other_collisions;

        debug_assert!(other_collisions.is_empty() ||
                      other_collisions[0].subject > src_ix as u16);

        let mut spawn = SpawnList::new();
        ExtendedObject::from(&object).update(
            self, &mut spawn, this_collisions,
            &self.writer.prev_snapshot);

        if unlikely!(1 != spawn.len()) {
            self.ensure_relative_object_capacity(spawn.len(), 1);
        }

        for object in spawn {
            self.write_object(object);
        }
    }

    fn ensure_relative_object_capacity(&mut self, add: usize, subtract: usize) {
        let new_count = self.reserved + add - subtract;
        if new_count > self.capacity {
            // Need to allocate more space. Target double what we immediately
            // need.
            let new_cap = new_count * 2;
            let mut new_objects = Vec::with_capacity(new_cap);
            let mut new_rtree = Vec::with_capacity(new_cap + 2);
            // Unsafe slicing here is a hard requirement since the length on
            // the input vectors is not set.
            new_objects.extend_from_slice(unsafe {
                self.objects.0.get_unchecked(..self.reserved)
            });
            new_rtree.extend_from_slice(unsafe {
                self.objects.1.get_unchecked(..2+self.reserved)
            });

            self.reserved = new_count;
            self.capacity = new_cap;
            self.objects = Arc::new((new_objects, new_rtree));
            self.writer.objects.set(self.objects.clone());
        }
    }

    /// Spontaneously spawn the given object.
    ///
    /// This is not used by the normal update pipeline as it is less efficient,
    /// but rather by external code to add objects to the snapshot.
    ///
    /// It is not public since it is only to be invoked in very particular
    /// circumstances.
    fn add_object(&mut self, object: CommonObject) {
        self.ensure_relative_object_capacity(1, 0);
        self.write_object(object);
    }

    #[inline]
    fn write_object(&mut self, object: CommonObject) {
        let ix = self.result.num_objects;
        debug_assert!(ix < self.capacity);

        let subtree_size = self.tree_builder.curr_size() as usize;
        self.tree_builder.push_next();
        // Get the left child of the next object into the cache. (Note that we
        // already advanced the builder past `ix`.) The `- 1` part is not
        // present since `ix` has not been incremented.
        unsafe {
            prefetch_read_data(self.objects.1.get_unchecked(
                2 + ix - self.tree_builder.curr_size() as usize), 3);
        }

        // Write the updated data into the new snapshot
        unsafe {
            ptr::write(self.objects.0.as_ptr().offset(ix as isize)
                       as *mut CommonObject, object);

            // Load the left child of the current bounding rhombus now,
            // to give the optimiser freedom to make the rest of the
            // code branch-free or lower it into the branch depending
            // on how efficient the union code is on the architecture.
            let left_rhombus_branch = *self.objects.1.get_unchecked(
                2 + ix - subtree_size - 1);
            let mut rhombus = object.bounding_rhombus();
            if 0 != subtree_size {
                rhombus = rhombus.union(self.prev_bounding_rhombus)
                    .union(left_rhombus_branch);
            }
            self.prev_bounding_rhombus = rhombus;

            ptr::write(self.objects.1.as_ptr().offset(2 + ix as isize)
                       as *mut BoundingRhombus, rhombus);
        }

        self.result.num_objects += 1;
    }

    /// Allocate data in the new snapshot's arena.
    ///
    /// Allocations are write-once. `data` is used to indicate both the size of
    /// the allocation as well as the actual data to write.
    pub fn alloc(&mut self, data: &[i32x4]) -> u16 {
        assert!(data.len() + self.result.arena_size <= 65536);

        let base = self.result.arena_size;
        self.result.arena_size += data.len();

        if self.result.arena_size > self.arena.capacity {
            let new_arena = ArenaMemory::with_capacity(
                2 * self.result.arena_size);
            unsafe {
                ptr::copy_nonoverlapping(self.arena.base, new_arena.base, base);
            }
            self.arena = Arc::new(new_arena);
            self.writer.arena.set(self.arena.clone());
        }

        unsafe {
            ptr::copy_nonoverlapping(
                data.as_ptr(), self.arena.base.offset(base as isize),
                data.len());
        }

        base as u16
    }

    /// Allocate a new object id.
    pub fn alloc_id(&mut self) -> u16 {
        self.result.free_ids.alloc()
    }

    /// Poll for events belonging to the given object.
    ///
    /// Events are returned in canonical order.
    pub fn poll_events<'b>(&'b self, id: u16)
                           -> impl 'b + Iterator<Item = &'b ObjectEvent> {
        lazy_static! {
            static ref EMPTY: BTreeSet<ObjectEvent> =
                BTreeSet::new();
        }

        self.events.get(&id).unwrap_or_else(|| &*EMPTY).iter()
    }
}

impl SnapshotWriter {
    /// Fully run the collision stage of the pipeline.
    ///
    /// This must be run after or concurrently with `run_update` or else it
    /// will deadlock.
    ///
    /// Returns the partial result of this pipeline stage.
    #[inline(never)]
    pub fn run_collision(&self) -> SnapshotCollisionResult {
        assert!(!self.update_started.swap(true, Ordering::Relaxed));

        let mut result = SnapshotCollisionResult {
            collisions: Vec::new(),
        };

        // Start from index 1 since object 0 has nothing to collide with
        let mut consumed = 1;
        // The scanner always points to the object *behind* consumed.
        let mut scanner = ImplicitTreeScanner::new();

        loop {
            // Spin-wait for more objects to be available
            let available = self.num_objects.load(Ordering::Acquire);
            if unlikely!(available <= consumed) {
                if self.update_done.load(Ordering::Acquire) {
                    break;
                }

                hint_core_should_pause();
                continue;
            }

            run_collision_through_available(
                &mut result, &mut consumed, &mut scanner,
                available, self.objects.get(), self.arena.get());
        }

        result
    }

}

fn run_collision_through_available(
    result: &mut SnapshotCollisionResult,
    consumed: &mut usize, scanner: &mut ImplicitTreeScanner,
    available: usize,
    objects: Arc<(Vec<CommonObject>, Vec<BoundingRhombus>)>,
    arena: Arc<ArenaMemory>
) {
    let objects_ptr = unsafe { slice_unchecked(&objects.0, 0, available) };
    while *consumed < available {
        let self_obj = *unsafe { get_unchecked(objects_ptr, *consumed) };
        let self_inv_rotate = Affine2dH::rotate_hex(-self_obj.theta());
        let bound = self_obj.bounding_rhombus().invert();

        // Search objects before this one for possible collisions.
        for (ix, candidate) in find_intersecting(
            scanner.navigate(),
            unsafe {
                // The slice excludes this object, which ensures we
                // don't waste time colliding with self.
                slice_unchecked(&objects.0, 0, *consumed)
            },
            unsafe {
                slice_unchecked(&objects.1, 0, *consumed)
            },
            bound
        ) {
            // Ignore same collision group
            if 0 != self_obj.collision_group() &&
                self_obj.collision_group() ==
                candidate.collision_group()
            { continue; }

            // Precise collision test
            let collision = check_precise_collision(
                self_obj, candidate, self_inv_rotate, &arena);

            if unlikely!(collision) {
                // Add the collision between this pair
                result.collisions.push(ImpendingCollision {
                    subject: *consumed as u16,
                    other: ix as u16,
                });
                result.collisions.push(ImpendingCollision {
                    subject: ix as u16,
                    other: *consumed as u16,
                });
            }
        }

        *consumed += 1;
        scanner.next();
    }
}

fn check_precise_collision(
    self_obj: CommonObject, candidate: CommonObject,
    self_inv_rotate: Affine2dH,
    arena: &ArenaMemory,
) -> bool {
    match (0 == self_obj.rounded_span(), 0 == candidate.rounded_span()) {
        // Point particle + point particle = no collision
        (true, true) => false,
        (false, true) => unsafe {
            let self_composite = CompositeObject::wrap(
                arena.decompress_ptr(
                    self_obj.extended_data(),
                    self_obj.data_dst_size()));
            let rel_pos = self_inv_rotate *
                (candidate.pos() - self_obj.pos()).dual();
            self_composite.test_point_collision(rel_pos.single()).is_some()
        },
        (true, false) => unsafe {
            let cand_composite = CompositeObject::wrap(
                arena.decompress_ptr(
                    candidate.extended_data(),
                    candidate.data_dst_size()));
            let rel_pos = Affine2dH::rotate_hex(-candidate.theta()) *
                (self_obj.pos() - candidate.pos()).dual();
            cand_composite.test_point_collision(rel_pos.single()).is_some()
        },
        (false, false) => unsafe {
            let self_composite = CompositeObject::wrap(
                arena.decompress_ptr(
                    self_obj.extended_data(),
                    self_obj.data_dst_size()));
            let cand_composite = CompositeObject::wrap(
                arena.decompress_ptr(
                    candidate.extended_data(),
                    candidate.data_dst_size()));
            let mut set = CollisionSet::new();
            self_composite.test_composite_collision(
                &mut set, self_obj, &cand_composite, candidate,
                self_inv_rotate);
            !set.is_empty()
        },
    }
}

bitflags! {
    /// Flags passed to `SnapshotWriter::finalise()`.
    pub struct FinaliseFlags : u32 {
        /// Garbage-collect the arena into a new allocation.
        const GC        = 1 << 0;
        /// Re-sort the objects. This has material effect on the simulation,
        /// and so must be done on the same frames for all members.
        const SORT      = 1 << 1;
    }
}

impl SnapshotWriter {
    /// Finalise this snapshot writer using the results of the pipeline stages.
    ///
    /// Returns a new writer that can be used for the next frame. The snapshot
    /// for this frame can be obtained from `prev_snapshot` on the returned
    /// writer.
    pub fn finalise(self, update: SnapshotUpdateResult,
                    mut collision: SnapshotCollisionResult,
                    flags: FinaliseFlags)
                    -> SnapshotWriter {
        let num_objects = update.num_objects;
        let objects = self.objects.get();
        drop(self.objects);

        let (mut objects, mut rhombus_tree) = match Arc::try_unwrap(objects) {
            Ok(v) => v,
            Err(_) => panic!("Leaked reference to SnapshotWriter::objects"),
        };

        unsafe {
            objects.set_len(num_objects);
            rhombus_tree.set_len(num_objects + 2);
        }

        let mut arena = self.arena.get();
        let mut arena_size = update.arena_size;
        let mut free_ids = update.free_ids;

        if flags.contains(FinaliseFlags::GC) {
            let mut free_ids_builder = IdFreeListBuilder::new();
            let new_arena = ArenaMemory::with_capacity(
                min(update.arena_size * 5 / 4, 65536));
            let mut new_arena_size = 0;
            assert!(new_arena.capacity >= update.arena_size);

            for object in &mut objects {
                if object.has_extended_data() {
                    let new_base = new_arena_size;
                    let old_base = object.extended_data() as usize;
                    let len = object.data_dst_size() as usize;

                    unsafe {
                        ptr::copy_nonoverlapping(
                            arena.base.offset(old_base as isize),
                            new_arena.base.offset(new_base as isize),
                            len);
                    }
                    new_arena_size += len;

                    object.set_extended_data(new_base as u16);
                }

                free_ids_builder.see(object.id());
            }

            arena = Arc::new(new_arena);
            arena_size = new_arena_size;
            free_ids = free_ids_builder.build();
        }

        if flags.contains(FinaliseFlags::SORT) {
            // u64s with the index of the object in the lower 16 bits and the
            // Hilbert coordinate in the upper 48.
            let mut refs = Vec::with_capacity(objects.len());
            refs.extend(objects.iter().enumerate().map(|(ix, obj)| {
                let mut hilbert_coord = xy_to_hilbert(obj.pos().repr());
                // High-speed objects can cause the tree quality to rapidly
                // deteriorate during non-sorted frames since such objects will
                // end up far away from objects they should be "near by". To
                // alleviate this somewhat, put "fast" objects at the end of
                // the array in no particular order. This does render the
                // top-levels of the tree mostly useless since the objects are
                // even farther away from their neighbours, but the only
                // objects that need to traverse this part of the tree during
                // collision detection are the high-speed objects themselves.
                //
                // Consider anything with a speed greater than one screen per
                // 5 screens/sec along both axes combined (L1 norm) to be
                // "fast".
                if obj.vx4().repr().abs().hsum_2() >=
                    4 * 5 * SCREEN * SECOND as i32
                {
                    hilbert_coord = !0;
                }

                (hilbert_coord << 16) | ix as u64
            }));
            // All values are unique, so we can use unstable sort if it's
            // faster.
            refs.sort_unstable();

            // Permute the object array accordingly
            let mut new_objs = Vec::with_capacity(objects.len());
            new_objs.extend(refs.iter().map(|&ix| {
                unsafe { get_unchecked(&objects, ix as usize & 0xFFFF) }
            }));
            objects = new_objs;

            // Compute the new collision tree
            let mut builder = ImplicitTreeBuilder::new();
            for (ix, obj) in objects.iter().enumerate() {
                let mut rhombus = obj.bounding_rhombus();
                let size = builder.curr_size() as usize;

                if size > 0 {
                    rhombus = rhombus
                        .union(*unsafe {
                            get_unchecked(&rhombus_tree, 2 + ix - size)
                        })
                        .union(*unsafe {
                            get_unchecked(&rhombus_tree, 2 + ix - 1)
                        });
                }

                unsafe {
                    *get_unchecked_mut(&mut rhombus_tree, ix + 2) = rhombus;
                }
                builder.push_next();
            }

            // Need to rewrite the collision pairs.
            // First, map from old index to new index.
            let mut indices = vec![0u16; refs.len()];
            for (ix, &r) in refs.iter().enumerate() {
                unsafe {
                    *get_unchecked_mut(&mut indices, r as usize & 0xFFFF) =
                        ix as u16;
                }
            }
            for c in &mut collision.collisions {
                c.subject = indices[c.subject as usize];
                c.other = indices[c.other as usize];
            }
        }

        collision.collisions.sort_unstable();

        let snapshot = Snapshot {
            objects, rhombus_tree,
            rhombus_tree_navigator:
                ImplicitTreeNavigator::root_of(num_objects as u32),
            arena: arena,
            arena_size: arena_size,
            collisions: collision.collisions,
            free_ids: free_ids,
            tick: self.prev_snapshot.tick + 1,
        };

        SnapshotWriter::wrap(Arc::new(snapshot), num_objects + 256)
    }
}
