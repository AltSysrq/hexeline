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

//! Definitions and functions for working with composite objects.
//!
//! All composite objects share several properties:
//!
//! - `CommonObject::extended_data()` is a pointer into the arena.
//!
//! - `CommonObject::data_dst_size()` gives the length of the array behind that
//!   pointer.
//!
//! - The allocated data is prefixed with a `CompositeObject`.
//!
//! The nominal (A,B) coordinate of a composite is always its centre of mass,
//! rather than the coordinate of the cell addressed as (0,0).
//!
//! The common composite data provides a few additional physical properties as
//! well as a bitset of which cells are populated. The bitset is addressed by
//! the integer (A,B) coordinates of each cell, with A being used as a row
//! index and B as a column index.

use std::borrow::Borrow;
use std::cmp::max;
use std::fmt;
use std::i16;
use std::ptr;

use simd::*;
use smallvec::{Array, SmallVec};

use physics::coords::*;

/// Unpacked form of `CompositeHeader`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct UnpackedCompositeHeader {
    /// The A coordinate of the (0,0) cell when the object itself is at (0,0)
    /// with theta=0.
    pub a_offset: i32,
    /// The B coordinate of the (0,0) cell when the object itself is at (0,0)
    /// with theta=0.
    pub b_offset: i32,
    /// Amount to subtract from the A coordinate of a cell to determine its row
    /// index. In other words, the minimum A coordinate of any row.
    pub row_offset: i16,
    /// The power of two of the number of rows in this composite.
    pub rows: u8,
    /// The power of two of the width, in cells, of each row.
    pub pitch: u8,
    /// The mass of this object.
    pub mass: u16,
}

impl UnpackedCompositeHeader {
    pub fn pack(self) -> CompositeHeader {
        CompositeHeader(i32x4::splat(0))
            .with_a_offset(self.a_offset)
            .with_b_offset(self.b_offset)
            .with_row_offset(self.row_offset)
            .with_rows(self.rows)
            .with_pitch(self.pitch)
            .with_mass(self.mass)
    }
}

/// The content of the first word in `CompositeObject`.
///
/// It specifies the size of the object as well as additional physical
/// properties.
#[derive(Clone, Copy)]
pub struct CompositeHeader(pub i32x4);

impl CompositeHeader {
    packed_field!(0:0[ 0..31]: i32 a_offset, set_a_offset, with_a_offset);
    packed_field!(0:1[ 0..31]: i32 b_offset, set_b_offset, with_b_offset);
    packed_field!(0:2[ 0..15]: i16 row_offset, set_row_offset, with_row_offset);
    packed_field!(0:2[16..23]:  u8 rows, set_rows, with_rows);
    packed_field!(0:2[24..31]:  u8 pitch, set_pitch, with_pitch);
    packed_field!(0:3[ 0..15]: u16 mass, set_mass, with_mass);

    /// Return the (a_offset, b_offset) vector.
    #[inline(always)]
    pub fn offset(self) -> Vhs {
        Vhs::from_repr(self.0)
    }

    pub fn unpack(self) -> UnpackedCompositeHeader {
        UnpackedCompositeHeader {
            a_offset: self.a_offset(),
            b_offset: self.b_offset(),
            row_offset: self.row_offset(),
            rows: self.rows(),
            pitch: self.pitch(),
            mass: self.mass(),
        }
    }
}

impl fmt::Debug for CompositeHeader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(&self.unpack(), f)
    }
}

/// The common prefix data shared by all composite objects.
///
/// The first word in a composite is always a `CompositeHeader`.
///
/// This is followed by the population bitset. A given cell (A,B) is found at
/// bit `((A & ((1 << rows) - 1)) << pitch) + (B & ((1 << pitch) - 1))`. The
/// use of modular arithmetic means that no bounds checks are necessary in
/// exchange for false positives, allowing the bounds checks to be moved until
/// after testing the bitset.
///
/// This is followed by the column offset array, which has a length equal to
/// `rows`. If the pitch is less than 128 (i.e., `pitch < 7`), each value is an
/// `i8`, packed 16 per word; otherwise, each value is an `i16`, packed 8 per
/// word. Each value in the column offset array gives the minimum B coordinate
/// of any populated cell in that row; the valid range for each row is thus the
/// `column_offset` to `column_offset + (1 << pitch) - 1`.
///
/// The bitset is padded to a multiple of 16 bits.
#[derive(Clone, Copy)]
pub struct CompositeObject<T : Borrow<[i32x4]>>(T);

impl<T : Borrow<[i32x4]>> CompositeObject<T> {
    /// Returns the header of this composite.
    #[inline(always)]
    pub fn header(&self) -> CompositeHeader {
        CompositeHeader(*unsafe { self.0.borrow().get_unchecked(0) })
    }

    #[inline(always)]
    fn bitset_start(&self) -> *const u8 {
        &self.0.borrow()[1] as *const i32x4 as *const u8
    }

    #[inline(always)]
    fn col_offset_array_start(&self) -> *const () {
        let header = self.header();
        (unsafe {
            self.bitset_start().offset(
                ((1 << header.rows() << header.pitch()) + 15) / 16 * 2)
        }) as *const ()
    }

    #[inline(always)]
    fn use_16bit_col_offset(&self) -> bool {
        self.header().pitch() >= 7
    }

    /// Return whether the cell at `(a,b)` is populated.
    ///
    /// If `(a,b)` is not in range, the result is unspecified but safe.
    #[inline(always)]
    pub fn is_populated(&self, a: i16, b: i16) -> bool {
        let header = self.header();
        let row = (a & ((1 << header.rows()) - 1)) as usize;
        let col = (b & ((1 << header.pitch()) - 1)) as usize;
        let bit = (row << header.pitch()) + col;
        let byte = bit >> 3;
        let shift = bit & 7;

        1 == (unsafe {
            ptr::read(self.bitset_start().offset(byte as isize))
        } >> shift) & 1
    }

    /// Return whether `a` is an in-bounds cell row for this composite.
    #[inline(always)]
    pub fn is_in_a_bound(&self, a: i16) -> bool {
        let header = self.header();
        let row = (a - header.row_offset()) as usize;
        row < 1 << header.rows()
    }

    /// Return whether `(a,b)` is an in-bounds cell address for this composite.
    #[inline(always)]
    pub fn is_in_bounds(&self, a: i16, b: i16) -> bool {
        let header = self.header();
        let row = (a - header.row_offset()) as usize;
        if row >= (1 << header.rows()) { return false; }

        let array = self.col_offset_array_start();
        let col_offset = if self.use_16bit_col_offset() {
            unsafe {
                ptr::read((array as *const i16).offset(row as isize))
            }
        } else {
            (unsafe {
                ptr::read((array as *const i8).offset(row as isize))
            }) as i16
        };

        let col = (b - col_offset) as usize;
        col < (1 << header.pitch())
    }

    /// Test for a collision with a point particle.
    ///
    /// `relative_pos` is the relative position of the point particle relative
    /// to this object's nominal coordinates, including adjusting for rotation.
    ///
    /// If there is a collision, returns the index of the cell affected.
    #[inline(always)]
    pub fn test_point_collision(&self, relative_pos: Vhs) -> Option<(i16,i16)> {
        let grid_relative_pos = relative_pos - self.header().offset();
        let (ia, ib) = grid_relative_pos.to_index();

        if self.is_in_a_bound(ia as i16) &&
            self.is_populated(ia as i16, ib as i16) &&
            self.is_in_bounds(ia as i16, ib as i16) &&
            ia >= i16::MIN as i32 && ia <= i16::MAX as i32 &&
            ib >= i16::MIN as i32 && ia <= i16::MAX as i32
        {
            Some((ia as i16, ib as i16))
        } else {
            None
        }
    }
}

impl<A : Array<Item = i32x4>> CompositeObject<SmallVec<A>> {
    /// Create a `CompositeObject` with the given header and cell bitmap.
    ///
    /// `cells` is called to yield an iterator over the `(A,B)` coordinates of
    /// populated cells. It must yield cells sorted ascending and without
    /// duplicates, and must always produce the same results.
    ///
    /// `pitch`, `rows`, and `row_offset` of `base` are adjusted automatically.
    ///
    /// ## Unsafety
    ///
    /// Behaviour is undefined if `cells` violates its general contract.
    pub unsafe fn build<T : Iterator<Item = (i16,i16)>, F : Fn () -> T>
        (mut base: CompositeHeader, cells: F) -> Self
    {
        let min_row = cells().next().expect(
            "Attempted to make composite with no cells").0;
        let mut max_row = min_row;
        let mut prev = (i16::MIN, i16::MIN);
        let mut row_start = i16::MIN;
        let mut max_row_width = 0;
        let mut max_abs_b = 0;

        for (a, b) in cells() {
            debug_assert!((a, b) > prev);

            if a != prev.0 {
                max_row_width = max(max_row_width, 1 + prev.1 - row_start);
                row_start = b;
            }

            max_abs_b = max(max_abs_b, b.abs());
            prev = (a, b);
            max_row = a;
        }
        max_row_width = max(max_row_width, 1 + prev.1 - row_start);

        fn log2_up(v: i16) -> u8 {
            (16 - (v - 1).leading_zeros()) as u8
        }

        let pitch = max(log2_up(max_row_width),
                        if max_abs_b <= 127 { 0 } else { 7 });
        let rows = log2_up(1 + max_row - min_row);
        let column_offset_fmt = if pitch >= 7 { 1 } else { 0 };
        base.set_row_offset(min_row);
        base.set_pitch(pitch);
        base.set_rows(rows);

        let mut dst = SmallVec::<A>::new();
        let data_bytes = ((1 << pitch << rows) + 15) / 16 * 2 +
            (1 << rows << column_offset_fmt);
        let capacity = 1 + (data_bytes + 15) / 16;
        dst.reserve(capacity);
        dst.push(base.0);
        while dst.len() < capacity {
            dst.push(i32x4::splat(0));
        }

        let bitset_base = &mut dst[1] as *mut i32x4 as *mut u8;
        let mut column_off = (/*unsafe*/ {
            bitset_base.offset(((1 << pitch << rows) + 15) / 16 * 2)
        }) as *mut ();

        prev = (min_row - 1, i16::MIN);
        for (a, b) in cells() {
            debug_assert!((a, b) > prev);
            // Need a loop rather than `if` to deal with gaps. The column
            // offsets we write for gaps are arbitrary so just use whatever is
            // convenient.
            while prev.0 < a {
                if 0 == column_offset_fmt {
                    /*unsafe*/ {
                        ptr::write(column_off as *mut i8, b as i8);
                        column_off = (column_off as *mut i8)
                            .offset(1) as *mut ();
                    }
                } else {
                    /*unsafe*/ {
                        ptr::write(column_off as *mut i16, b);
                        column_off = (column_off as *mut i16)
                            .offset(1) as *mut ();
                    }
                }
                prev.0 += 1;
            }
            prev = (a, b);

            let bit_row = (a & ((1 << rows) - 1)) as usize;
            let bit_col = (b & ((1 << pitch) - 1)) as usize;
            let bit = (bit_row << pitch) + bit_col;
            let byte = bit >> 3;
            let shift = bit & 7;

            /*unsafe*/ {
                let bptr = bitset_base.offset(byte as isize);
                ptr::write(bptr, (1 << shift) | ptr::read(bptr));
            }
        }

        // Since the rows count is rounded up, there may be rows whose columns
        // offsets weren't set. That's fine; with no cells, any value is fine,
        // and we zero-initialised everything already.

        debug_assert!((column_off as usize) <=
                      &dst[dst.len() - 1] as *const i32x4 as usize + 16,
                      "Wrote too many rows (got to {:?}, \
                       allocation ends at {:?}+16", column_off,
                      &dst[dst.len() - 1] as *const i32x4);

        CompositeObject(dst)
    }
}

#[cfg(test)]
mod test {
    use proptest;
    use proptest::strategy::Strategy;

    use super::*;

    proptest! {
        #[test]
        fn packed_header_preserves_all_fields(
            a_offset in proptest::num::i32::ANY,
            b_offset in proptest::num::i32::ANY,
            row_offset in proptest::num::i16::ANY,
            rows in proptest::num::u8::ANY,
            pitch in proptest::num::u8::ANY,
            mass in proptest::num::u16::ANY
        ) {
            let orig = UnpackedCompositeHeader {
                a_offset, b_offset, row_offset,
                rows, pitch, mass,
            };
            assert_eq!(orig.pack().unpack(), orig);
        }

        #[test]
        fn builds_correct_bitset(
            ref cells in proptest::collection::btree_set(
                (-64i16..64i16, -64i16..64i16), 1..64)
            .prop_union(
                proptest::collection::btree_set(
                    (-1023i16..1024i16, -1023i16..1024i16), 1..1024))
        ) {
            let object = unsafe {
                CompositeObject::<SmallVec<[i32x4;32]>>::build(
                    UnpackedCompositeHeader::default().pack(),
                    || cells.iter().map(|&v| v))
            };

            // All cells are marked populated
            for &(a, b) in cells {
                assert!(object.is_populated(a, b),
                        "Cell ({},{}) not marked populated", a, b);
                assert!(object.is_in_bounds(a, b),
                        "Cell ({},{}) is out of bounds", a, b);
            }

            // Coordinates not covered by `cells` are not populated or are out
            // of bounds
            for a in -128i16..128i16 {
                for b in -128i16..128i16 {
                    if !cells.contains(&(a, b)) {
                        let pop = object.is_populated(a, b);
                        let in_bounds = object.is_in_bounds(a, b);
                        assert!(!pop || !in_bounds,
                                "Non-cell ({},{}) marked present; \
                                 populated = {}, in_bounds = {}",
                                a, b, pop, in_bounds);
                    }
                }
            }
        }
    }
}
