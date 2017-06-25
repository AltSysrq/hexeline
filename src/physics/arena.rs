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

/*!
Facilities for allocating and managing extended object data.

`CommonObject` has an `extended_data` value which is a compressed pointer into
an arena. Object types which need more than 24 bits of extended data must
allocate space in the arena to store this information.

An arena is a single block of contiguous memory, addressed by indexing `i32x4`
values. Given 16-bit compressed pointers, we thus have a maximum of 1MB of
space. Memory is allocated in a strictly linear fashion, using a simple
incrementing allocator. Memory below the "eden" address is immutable; eden
advances each time the arena is frozen. This means multiple snapshots can share
the same arena.

However, it also means that most mutations need to allocate new memory, and
there is no way to free memory other than freeing the whole arena. Therefore,
when the arena starts getting full, it must be evacuated to a fresh arena by
what is essentially a moving garbage collection process.

Note that a compressed pointer of 0 is perfectly valid. However, since
allocations are created in strictly ascending order, in contexts where there is
guaranteed to be at least one allocation before a particular pointer (e.g.,
when a pointer is itself stored in such an allocation), the value 0 is still
useful as a "null pointer".
*/

use std::mem;
use std::ops;
use std::slice;
use std::sync::Arc;

use simd::i32x4;

/// Manages the raw memory allocation behind an arena.
///
/// This is essentially a `Vec` in disguise.
///
/// We use a dynamic memory expansion strategy by at least doubling `capacity`
/// when space is exhausted. This allows us to avoid allocating the full 1MB in
/// cases where it can be avoided. Note that these expansions do not entail an
/// evacuation, and end up bringing obsolete data long, since the arena has no
/// way to determine what data is live, much less update pointers to it.
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

/// A handle on an arena which the holder cannot mutate.
///
/// This is "const" in that it is immutable as long as it is held, but will not
/// necessarily be immutable indefinitely, and so it is not `Clone`.
#[derive(Debug)]
pub struct ConstArena {
    mem: Arc<ArenaMemory>,
    base: *mut i32x4,
    len: usize,
}

impl ConstArena {
    fn clone_internal(&self) -> Self {
        ConstArena {
            mem: self.mem.clone(),
            base: self.base,
            len: self.len,
        }
    }
}

/// An arena which has been fully frozen.
///
/// The contents are guaranteed to be immutable indefinitely, and so it can be
/// freely cloned or sent between threads.
#[derive(Debug)]
pub struct FrozenArena(ConstArena);
unsafe impl Send for FrozenArena { }
unsafe impl Sync for FrozenArena { }

impl ops::Deref for FrozenArena {
    type Target = ConstArena;

    fn deref(&self) -> &ConstArena {
        &self.0
    }
}

impl Clone for FrozenArena {
    fn clone(&self) -> Self {
        FrozenArena(self.0.clone_internal())
    }
}

/// An arena which is writable.
///
/// There are two ways to write to the arena:
///
/// - Allocate new memory and write to it.
///
/// - Overwrite memory above the "eden" mark.
///
/// The eden mark is the greatest address that any `FrozenArena` may include.
/// Memory below the eden mark cannot be mutated; memory above it can be, since
/// only this `ActiveArena` can see it.
#[derive(Debug)]
pub struct ActiveArena {
    base: ConstArena,
    eden: usize,
    limit: usize,
}

impl ops::Deref for ActiveArena {
    type Target = ConstArena;

    fn deref(&self) -> &ConstArena {
        &self.base
    }
}

impl ConstArena {
    /// Decompress the given pointer to access the memory behind it.
    ///
    /// ## Unsafety
    ///
    /// No bounds checks are performed if debugging assertions are disabled.
    #[inline]
    pub unsafe fn decompress(&self, offset: u16, len: u8) -> &[i32x4] {
        debug_assert!((offset as usize) + (len as usize) < self.len);
        slice::from_raw_parts(self.base.offset(offset as isize),
                              len as usize)
    }

    /// Returns an upper bound on the number of `i32x4` values reachable from
    /// this `ConstArena`.
    pub fn used(&self) -> usize {
        self.len
    }

    /// Returns the number of `i32x4` values in the array backing this
    /// `ConstArena`.
    pub fn capacity(&self) -> usize {
        self.mem.capacity
    }
}

impl ActiveArena {
    /// Create a new, empty arena with space for `capacity` elements.
    pub fn with_capacity(capacity: usize) -> Self {
        let mem = ArenaMemory::with_capacity(capacity);
        let base = mem.base;
        ActiveArena {
            base: ConstArena {
                mem: Arc::new(mem),
                base: base,
                len: 0,
            },
            eden: 0,
            limit: capacity,
        }
    }

    /// Return a new, empty arena with a capacity determined from the current
    /// usage of this arena.
    pub fn fresh(&self) -> Self {
        let capacity = if self.base.len <= self.limit / 2 {
            self.limit / 2
        } else {
            self.limit
        };
        ActiveArena::with_capacity(capacity)
    }

    /// Like `decompress()`, but returns mutable memory.
    ///
    /// ## Unsafety
    ///
    /// No bounds checks are performed if debugging assertions are disabled.
    #[inline]
    pub unsafe fn decompress_mut(&mut self, offset: u16, len: u8)
                                 -> &mut [i32x4] {
        debug_assert!(offset as usize >= self.eden);
        debug_assert!((offset as usize) + (len as usize) < self.base.len);
        slice::from_raw_parts_mut(self.base.base.offset(offset as isize),
                                  len as usize)
    }

    /// Allocate space for `data`, copy `data` to that space, and return a
    /// compressed pointer referring to it.
    ///
    /// ## Panics
    ///
    /// Panics if allocating `data` would grow the arena to more than 65536
    /// elements.
    pub fn alloc(&mut self, data: &[i32x4]) -> u16 {
        if self.base.len + data.len() > self.limit {
            let min_len = self.base.len + data.len();
            self.grow(min_len);
        }

        let offset = self.base.len;
        self.base.len += data.len();
        unsafe {
            self.decompress_mut(offset as u16, data.len() as u8)
                .copy_from_slice(data);
        }
        offset as u16
    }

    /// Return a compressed pointer above the eden mark pointing to the same
    /// data that `offset` does.
    ///
    /// If `offset` is above the eden mark, returns `offset`. Otherwise,
    /// allocates `len` elements, copies them from `offset`, and returns the
    /// new allocation.
    #[inline]
    pub fn realloc(&mut self, offset: u16, len: u8) -> u16 {
        if (offset as usize) >= self.eden {
            offset
        } else {
            let old = unsafe {
                &(*(self.decompress(offset, len) as *const [i32x4]))[..]
            };
            self.alloc(old)
        }
    }

    /// Freeze this arena.
    ///
    /// The current state of the arena is returned as a `FrozenArena`. A new
    /// `ActiveArena` is also returned which can be used to continue appending
    /// to the same memory.
    pub fn freeze(self) -> (FrozenArena, ActiveArena) {
        let new_eden = self.base.len;
        (FrozenArena(self.base.clone_internal()),
         ActiveArena {
            base: self.base, eden: new_eden,
            limit: self.limit
        })
    }

    fn grow(&mut self, min_size: usize) {
        use std::cmp::{max, min};

        if min_size > 65536 {
            panic!("Tried to grow Arena beyond 65536 elements");
        }

        let new_size = min(65536, max(min_size, 2 * self.limit));
        self.base.mem = Arc::new(ArenaMemory::with_capacity(new_size));
        self.limit = new_size;
    }
}
