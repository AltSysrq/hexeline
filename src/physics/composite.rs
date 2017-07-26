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
use std::cmp::{max, min};
use std::fmt;
use std::i16;
use std::ptr;

use arrayvec::ArrayVec;
use simd::*;
use smallvec::{Array, SmallVec};

use intext::*;
use physics::common_object::*;
use physics::coords::*;
use physics::xform::Affine2dH;

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

pub const COLLISION_SET_SIZE: usize = 16;
pub type CollisionSet = ArrayVec<[[(i16,i16);2];COLLISION_SET_SIZE]>;

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
        let row = (a as i32 - header.row_offset() as i32) as usize;
        if row >= (1 << header.rows()) { return false; }

        let col_offset = unsafe { self.col_offset_unchecked(row) };

        let col = (b as i32 - col_offset as i32) as usize;
        col < (1 << header.pitch())
    }

    /// Returns an iterator over the logical row indices in this composite.
    pub fn rows(&self) -> impl Iterator<Item = i16> {
        let row_offset = self.header().row_offset() as i32;
        (0..(1i32 << self.header().rows())).map(
            move |row| (row + row_offset) as i16)
    }

    /// Returns the logical index of the first column in the given logical row.
    pub fn col_offset(&self, row: i16) -> i16 {
        let phys_row = (row as i32 - self.header().row_offset() as i32)
            as usize;
        assert!(phys_row <= 1 << self.header().rows());
        unsafe {
            self.col_offset_unchecked(phys_row)
        }
    }

    /// Returns an iterator over the logical indices of columns in the given
    /// logical row.
    pub fn cols_in_row(&self, row: i16) -> impl Iterator<Item = i16> {
        let col_offset = self.col_offset(row) as i32;
        (0..(1i32 << self.header().pitch())).map(
            move |col| ((col + col_offset) as i16))
    }

    /// Returns an iterator over all 2D cell indices considered "in-bounds" in
    /// this composite.
    pub fn indices<'a>(&'a self) -> impl 'a + Iterator<Item = (i16,i16)> {
        self.rows().flat_map(
            move |row| self.cols_in_row(row).map(move |col| (row, col)))
    }

    /// Returns an iterator over the logical indices of cells which are
    /// populated in this composite.
    pub fn cells<'a>(&'a self) -> impl 'a + Iterator<Item = (i16,i16)> {
        self.indices().filter(move |&(row, col)| self.is_populated(row, col))
    }

    #[inline(always)]
    unsafe fn col_offset_unchecked(&self, row: usize) -> i16 {
        let array = self.col_offset_array_start();
        if self.use_16bit_col_offset() {
            ptr::read((array as *const i16).offset(row as isize))
        } else {
            (ptr::read((array as *const i8).offset(row as isize))) as i16
        }
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
            (ia as i16 as i32) == ia && (ib as i16 as i32) == ib
        {
            Some((ia as i16, ib as i16))
        } else {
            None
        }
    }

    /// Computes an upper bound on the collision radius of this object (as per
    /// `CommonObject::rounded_radius`, but not yet rounded) from the centre of
    /// gravity.
    pub fn calc_radius(&self) -> u32 {
        // Sacrifice a bit of precision to avoid overflow
        // Dropping 4 bits means we can support a max distance of 2**(15+4) =
        // 2**19 (8 screens, far more than we'll ever need) rather than just
        // 2**15 (half a screen). We just need to adjust the final calculation
        // to be a bit more conservative.
        const SHIFT: u8 = 4;

        let header = self.header();
        let offset = header.offset() >> SHIFT;
        let mut max_dist = 0;

        // Search for the cell which is the farthest from the centre of
        // gravity. It could be in any row, so we must scan all the rows.
        for row in self.rows() {
            let mut col_iter = self.cols_in_row(row).filter(
                |&col| self.is_populated(row, col));
            if let Some(first_col) = col_iter.next() {
                let last_col = col_iter.last().unwrap_or(first_col);

                let first_pos = Vhs(row as i32, first_col as i32);
                let last_pos = Vhs(row as i32, last_col as i32);
                let first_pos = first_pos << (CELL_HEX_SHIFT - SHIFT);
                let last_pos = last_pos << (CELL_HEX_SHIFT - SHIFT);
                let first_pos = first_pos + offset;
                let last_pos = last_pos + offset;

                max_dist = max(first_pos.redundant().nsw_l2_squared(),
                               max(last_pos.redundant().nsw_l2_squared(),
                                   max_dist));
            }
        }

        // Compute actual distance, then round up the precision lost by
        // `SHIFT`, and add a cell radius to account for the cell coordinate
        // being the centre rather than the edge.
        CommonObject::radius_of_l2(
            ((max_dist.sqrt_up() + 1) << SHIFT) + CELL_L2_VERTEX as u32)
    }

    /// Compare two composites for collisions.
    ///
    /// `self_inverse_rot_xform` is the value of
    /// `Affine2dH::rotate_hex(-self_obj.theta())`.
    ///
    /// Returns the set of all pairs of colliding cells, or a subset of them if
    /// the vector is filled.
    ///
    /// Note that this operation is not commutative since the cell collision
    /// check is only approximate and due to precision limitations.
    pub fn test_composite_collision<R : Borrow<[i32x4]>>(
        &self, self_obj: CommonObject,
        that: &CompositeObject<R>, that_obj: CommonObject,
        self_inverse_rot_xform: Affine2dH
    ) -> CollisionSet {
        // To check the collision, we first take the bounding rhombus of `that`
        // and overlay it on the grid of `self` to find candidate cells; we
        // then check each cell of `self` within that rhombus against `that` to
        // see if anything overlaps.
        //
        // There is not actually any reason to try to compare the smaller
        // object against the larger one or vice-versa. If `self` is
        // substantially larger than `that`, the bounding rhombus of `that`
        // constrains us to a smaller number of cells. If `self` is
        // substantially smaller, even if `that` covers all of `self`, we still
        // don't have too many cells to check.

        // Determine where the centre of gravity of `that` falls within the
        // grid of `self`.
        //
        // Code below assumes this value stays in dual coordinates.
        let that_cog_self_grid: Vhd = self_inverse_rot_xform *
            (that_obj.pos() - self_obj.pos()).dual() -
            self.header().offset().dual();
        // Determine the bounding rhombus of `that` within the grid of `self`.
        // We can use a simple bitshift with a slop factor for a fast
        // approximation of what area of the grid is covered.
        let that_radius = (that_obj.rounded_radius() as i32)
            << ROUNDED_RADIUS_SHIFT;
        let that_bounds_self_grid = that_cog_self_grid.repr() +
            i32x4::new(-that_radius, -that_radius, that_radius, that_radius);
        let that_bounds_self_grid =
            (that_bounds_self_grid >> CELL_HEX_SHIFT) +
            i32x4::new(-1, -1, 1, 1);

        // Now determine the position of our (0,0) cell within `that`'s grid.
        let that_inverse_rot_xform = Affine2dH::rotate_hex(-that_obj.theta());
        let self_rot_xform = Affine2dH::rotate_hex(self_obj.theta());
        let origin_pos = self_obj.pos().dual() +
            self_rot_xform * self.header().offset().dual();
        let self_origin_that_grid =
            (that_inverse_rot_xform * (origin_pos - that_obj.pos().dual()))
            .single() - that.header().offset();

        // Determine the displacement within `that`'s grid for moving (+1,0)
        // and (0,+1) in our own cell grid
        let grid_displacement = that_inverse_rot_xform * self_rot_xform *
            Vhl(CELL_HEX_SIZE, 0, 0, CELL_HEX_SIZE);

        // Determine the first and last rows to scan, and start tracking the
        // base coordinate for that row.
        let first_row = max(that_bounds_self_grid.extract(0),
                            self.header().row_offset() as i32);
        let last_row = min(that_bounds_self_grid.extract(2),
                           (self.header().row_offset() as i32) +
                           (1 << self.header().rows()) - 1);
        let mut row_zero: Vhs = self_origin_that_grid +
            grid_displacement.fst() * Vhs(first_row, first_row);

        // Scan the columns of each row for overlapping cells
        let mut collisions = ArrayVec::new();
        let pitch = 1 << self.header().pitch();
        for row in first_row..last_row + 1 {
            let col_offset = unsafe {
                self.col_offset_unchecked(
                    (row - self.header().row_offset() as i32) as usize
                )
            };

            let first_col = max(that_bounds_self_grid.extract(1),
                                col_offset as i32);
            let last_col = min(that_bounds_self_grid.extract(3),
                               col_offset as i32 + pitch - 1);
            let mut coord: Vhs = row_zero + grid_displacement.snd() *
                Vhs(first_col, first_col);
            for col in first_col..last_col + 1 {
                if self.is_populated(row as i16, col as i16) {
                    let (rhs_a, rhs_b) = coord.to_grid_overlap();
                    macro_rules! check {
                        ($off:expr) => {
                            let a = rhs_a.extract($off);
                            let b = rhs_b.extract($off);
                            if that.is_populated(a as i16, b as i16) &&
                                that.is_in_bounds(a as i16, b as i16) &&
                                (a as i16 as i32) == a &&
                                (b as i16 as i32) == b
                            {
                                if collisions.push([(row as i16, col as i16),
                                                    (a as i16, b as i16)])
                                    .is_some()
                                {
                                    return collisions;
                                }
                            }
                        }
                    }
                    check!(0);
                    check!(1);
                    check!(2);
                    check!(3);
                }
                coord = coord + grid_displacement.snd();
            }

            row_zero = row_zero + grid_displacement.fst();
        }

        collisions
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

impl<T : Borrow<[i32x4]>> fmt::Debug for CompositeObject<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let header = self.header();

        let mut s = f.debug_struct("CompositeObject");
        s.field("header", &header);

        for row in self.rows() {
            if self.cols_in_row(row).any(|col| self.is_populated(row, col)) {
                s.field(&format!("row[{}]", row), &DebugRow(self, row));
            }
        }
        s.finish()
    }
}

struct DebugRow<T>(T, i16);
impl<'a, T : Borrow<[i32x4]>> fmt::Debug
for DebugRow<&'a CompositeObject<T>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let logical_row = self.1;
        let mut list = f.debug_list();
        for col in self.0.cols_in_row(logical_row) {
            if self.0.is_populated(logical_row, col) {
                list.entry(&col);
            }
        }
        list.finish()
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;
    use std::num::Wrapping;

    use proptest;
    use proptest::strategy::{BoxedStrategy, Strategy};

    use super::*;

    #[derive(Debug)]
    struct ArbComposite {
        common: CommonObject,
        data: CompositeObject<SmallVec<[i32x4;32]>>,
    }

    impl ArbComposite {
        fn cell_coord(&self, row: i16, col: i16) -> Vhd {
            let off = self.data.header().offset() +
                (Vhs(row as i32, col as i32) << CELL_HEX_SHIFT);
            self.common.pos().dual() +
                Affine2dH::rotate_hex(self.common.theta()) * off.dual()
        }
    }

    fn arb_composite() -> BoxedStrategy<ArbComposite> {
        (proptest::collection::btree_set(
            (-8i16..8i16, -8i16..8i16), 1..64),
         -65536..65536, -65536..65536,
         proptest::num::i16::ANY,
         -4096..4096, -4096..4096).prop_map(
            |(cells, a, b, theta, a_offset, b_offset)| {
                let mut ret = ArbComposite {
                    common: UnpackedCommonObject {
                        a, b,
                        theta: Wrapping(theta),
                        .. UnpackedCommonObject::default()
                    }.pack(),
                    data: unsafe { CompositeObject::build(
                        UnpackedCompositeHeader {
                            a_offset, b_offset,
                            .. UnpackedCompositeHeader::default()
                        }.pack(),
                        || cells.iter().map(|&v| v))
                    },
                };
                let rad = ret.data.calc_radius();
                ret.common.set_rounded_radius(
                    CommonObject::round_radius(rad));
                ret
            }).boxed()
    }

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

    proptest! {
        #![proptest_config(proptest::test_runner::Config {
            cases: 65536,
            .. proptest::test_runner::Config::default()
        })]

        #[test]
        fn composite_composite_collisions_roughly_correct(
            ref lhs in arb_composite(),
            ref rhs in arb_composite()
        ) {
            // For an exact solution, hexagons would collide somewhere between
            // 2*CELL_L2_EDGE and 2*CELL_L2_VERTEX. Increase the latter by
            // 13/10 to account for various precision loss. The former we
            // reduce to half of what it would normally be due to the way the
            // collision test is approximated.
            const AGGRESSIVE_DIST: u32 = CELL_L2_VERTEX as u32 * 2 * 13/10;
            const CONSERVATIVE_DIST: u32 = CELL_L2_EDGE as u32;

            let lhs_inv_rot = Affine2dH::rotate_hex(-lhs.common.theta());

            // Do a naïve n² collision check between all possible pairs.
            // Anything within `CONSERVATIVE_DIST` MUST be detected as a
            // collision; anything outside `AGGRESSIVE_DIST` MUST NOT be
            // detected as a collision. Things in between may or may not be
            // considered to collide.
            let mut aggressive_hits = HashSet::new();
            let mut conservative_hits = HashSet::new();
            for (lrow, lcol) in lhs.data.cells() {
                let lpos = lhs.cell_coord(lrow, lcol);
                for (rrow, rcol) in rhs.data.cells() {
                    let rpos = rhs.cell_coord(rrow, rcol);
                    let displacement = (lpos - rpos).redundant();
                    // Check LInf distance first since L2 distance could overflow
                    if displacement.linf() > AGGRESSIVE_DIST { continue; }

                    let l2d = displacement.nsw_l2_squared();
                    if l2d < CONSERVATIVE_DIST * CONSERVATIVE_DIST {
                        conservative_hits.insert([(lrow, lcol), (rrow, rcol)]);
                    }
                    if l2d <= AGGRESSIVE_DIST * AGGRESSIVE_DIST {
                        aggressive_hits.insert([(lrow, lcol), (rrow, rcol)]);
                    }
                }
            }

            let result = lhs.data.test_composite_collision(
                lhs.common, &rhs.data, rhs.common,
                lhs_inv_rot);
            for item in &result {
                assert!(aggressive_hits.contains(item),
                        "Detected {:?} as a collision, but they should not \
                         have collided", item);
            }
            // Only check that all conservative items were found if the result
            // array was not saturated
            if result.len() < COLLISION_SET_SIZE {
                for item in &conservative_hits {
                    assert!(result.iter().any(|i| i == item),
                            "Failed to detect {:?} as a collision", item);
                }
            }
        }
    }
}
