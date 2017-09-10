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

use std::sync::Arc;

use bit_set::BitSet;

/// Builder for `IdFreeList`.
#[derive(Debug, Clone)]
pub struct IdFreeListBuilder {
    available: BitSet,
}

/// Tracks what `u16` ids in a given space are free.
///
/// There is no way to explicitly mark an id as free. The id allocator simply
/// walks monotonically forward. Recovering unused ids requires a garbage
/// collection pass to build a new `IdFreeList`.
///
/// It is cheap to clone `IdFreeList` as they will share their underlying
/// memory.
#[derive(Debug, Clone)]
pub struct IdFreeList {
    available: Arc<BitSet>,
    next: usize,
}

impl IdFreeListBuilder {
    /// Initialise a new builder with all non-zero ids initially set as
    /// available.
    pub fn new() -> Self {
        let mut this = IdFreeListBuilder {
            available: BitSet::with_capacity(65536),
        };
        this.available.extend(1..65536);
        this
    }

    /// Marks the given id as in-use.
    pub fn see(&mut self, id: u16) {
        self.available.remove(id as usize);
    }

    /// Finalise this builder into an `IdFreeList`.
    pub fn build(self) -> IdFreeList {
        IdFreeList {
            available: Arc::new(self.available),
            next: 1,
        }
    }
}

impl IdFreeList {
    /// Allocate a new, unused id from the free list.
    pub fn alloc(&mut self) -> u16 {
        while self.available.contains(self.next) {
            self.next += 1;
        }

        debug_assert!(self.next < 65536);
        self.next as u16
    }

    /// Return a lower bound on the next id that will be generated.
    pub fn approx_next(&self) -> usize {
        self.next
    }
}
