//-
// Copyright (c) 2017, 2019 Jason Lingle
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

use std::borrow::{Borrow, BorrowMut};
use std::cmp::{max, min};
use std::fmt;
use std::i16;
use std::slice;

use arrayvec::ArrayVec;
use odds::slice_unchecked;
use simd::*;
use smallvec::{Array, SmallVec};

use simdext::*;
use intext::*;
use physics::common_object::*;
use physics::coords::*;
use physics::xform::{AFFINE_POINT, Affine2dH};

// Amount of precision to shift off of `AFFINE_POINT` in a couple places.
const SUBAFFINE_SHIFT: u32 = 4;

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
    /// The actual number of rows in this composite, excluding the empty first
    /// row and any padding rows after the last one.
    pub row_count: u16,
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
            .with_row_count(self.row_count)
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
    packed_field!(0:3[16..31]: u16 row_count, set_row_count, with_row_count);

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
            row_count: self.row_count(),
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

/// The number of cells packed into a single chunk.
///
/// Parts of the implementation depend strongly on this value, so don't change
/// it lightly.
pub const CHUNK_WIDTH: u32 = 24;

/// A chunk of the composite cell bitmap.
///
/// Each chunk represents up to `CHUNK_WIDTH` cells on the same row.
/// Additionally, it holds a bit for the `CHUNK_WIDTH+1`th cell on the row, as
/// well as the same range of cells on the next row (which may not be a
/// complete representation of that row).
///
/// Cell indices within a chunk are always even numbers between 0 and
/// CHUNK_WIDTH*2.
///
/// Additionally, each chunk stores the column offset for the row it is found
/// within.
#[derive(Clone, Copy)]
pub struct Chunk(i64);

impl Chunk {
    /// Returns the column base for this row.
    ///
    /// That is, the minimum B coordinate of any cell slot in this row.
    #[inline(always)]
    pub fn col_base(self) -> i16 {
        (self.0 >> 50) as i16
    }

    #[inline(always)]
    fn assert_valid_index(ix: u16) {
        debug_assert!(0 == ix & 1 && (ix as u32) <= CHUNK_WIDTH*2);
    }

    /// Returns the neighbourhood around the `ix`th cell in this chunk.
    #[inline(always)]
    pub fn neighbourhood(self, ix: u16) -> Neighbourhood {
        Chunk::assert_valid_index(ix);
        Neighbourhood((self.0 >> ix) as u32 & 0xF)
    }

    /// Returns whether the cell at the given index is populated.
    #[inline(always)]
    pub fn is_populated(self, ix: u16) -> bool {
        Chunk::assert_valid_index(ix);
        0 != (self.0 >> ix) & 1
    }

    /// Returns whether any cells belonging to this row are present in this
    /// chunk.
    #[inline(always)]
    pub fn has_any_on_row(self) -> bool {
        0 != self.0 & 0x555555555555
    }

    fn set_col0(&mut self, col: u16) {
        self.0 |= 1 << col;
    }

    fn set_col1(&mut self, col: u16) {
        self.0 |= 1 << col+1;
    }

    fn set_col_base(&mut self, base: i16) {
        self.0 |= (base as i64) << 50;
    }
}

impl fmt::Debug for Chunk {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Chunk({},{:048b})", self.col_base(),
               self.0 & ((1 << 2*CHUNK_WIDTH) - 1))
    }
}

/// A 2x2 neighbourhood of the cell bitmap.
#[derive(Clone, Copy, Debug)]
pub struct Neighbourhood(u32);

impl Neighbourhood {
    /// Returns whether any cells in the neighbourhood are present.
    #[inline(always)]
    pub fn any(self) -> bool {
        0 != self.0
    }

    /// Returns whether the cell at <0,0> is present.
    #[inline(always)]
    pub fn c00(self) -> bool {
        0 != self.0 & 1
    }

    /// Returns whether the cell at <1,0> is present.
    #[inline(always)]
    pub fn c10(self) -> bool {
        0 != self.0 & 2
    }

    /// Returns whether the cell at <0,1> is present.
    #[inline(always)]
    pub fn c01(self) -> bool {
        0 != self.0 & 4
    }

    /// Returns whether the cell at <1,1> is present.
    #[inline(always)]
    pub fn c11(self) -> bool {
        0 != self.0 & 8
    }

    /// Returns the neighbourhood of cells which are overlapped by the unit
    /// hexagon at `v`.
    pub fn hits(self, v: Vhs) -> Self {
        Neighbourhood(self.0 & v.to_grid_overlap_mask())
    }
}

/// A coordinate of a chunk and bit within a row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChunkRowCoord {
    /// The offset of the chunk within the row.
    pub chunk: u16,
    /// The bit (always a multiple of two) within the chunk.
    pub bit: u16,
}

fn prepare_rbi_chunk(chunk: Chunk) -> u64 {
    // Filter out everything other than the cells we care about, then put a
    // sentinel of 1 at the top (which will get right-shifted to 0 when the
    // scanner reaches the end).
    chunk.0 as u64 & 0x555555555555 | (1 << 2*CHUNK_WIDTH)
}

/// An iterator which returns offsets of present cells within a row.
///
/// The iterator is infinite; it simply wraps around every time it reaches the
/// end. This also means the last value is meaningless.
#[derive(Clone, Copy, Debug)]
struct RowBitIter<'a> {
    /// The data for this row. Always a power of two in length.
    row: &'a [Chunk],
    /// The index in `row` of the next chunk to load.
    next_chunk: usize,
    /// The current chunk. Prepared with `prepare_rbi_chunk`. This is
    /// essentially a bitset of all cells in this row at odd bit indices, with
    /// a 1 bit sentinel after the very last bit.
    current: u64,
}

impl<'a> RowBitIter<'a> {
    fn start_of_row<T : Borrow<[i32x4]>>(
        composite: &'a CompositeObject<T>,
        row: i16
    ) -> (Self, i16) {
        let row = composite.row(row);
        let first_chunk = unsafe { row.get_unchecked(0) };
        (Self::from(row, composite.col_index(first_chunk.col_base())),
         first_chunk.col_base())
    }

    /// Returns a partially initialised iterator positioned somewhere in the
    /// given row, and the index of the first column in that row.
    fn in_row<T : Borrow<[i32x4]>>(
        composite: &'a CompositeObject<T>,
        row: i16
    ) -> (Self, i16) {
        let row = composite.row(row);
        let ix = unsafe { row.get_unchecked(0) }.col_base();
        (RowBitIter { row, next_chunk: 0, current: 0 }, ix)
    }

    fn from(row: &'a [Chunk], col: ChunkRowCoord) -> Self {
        let mut this = RowBitIter {
            row,
            next_chunk: 0,
            current: 0,
        };
        this.reset(col);
        this
    }

    fn reset(&mut self, col: ChunkRowCoord) {
        self.next_chunk = (col.chunk as usize + 1) & (self.row.len() - 1);
        self.current = prepare_rbi_chunk(unsafe {
            *self.row.get_unchecked(col.chunk as usize)
        }) >> col.bit;
    }

    fn scan_through(self, base: i16, end: i16)
                    -> impl 'a + Iterator<Item = i16> {
        // Need to offset base by 1 since each iteration returns at least 1 but
        // if there's a cell in the first one, we want an output of `base`.
        self.scan(base - 1,
                  |accum, step| { *accum += step as i16; Some(*accum) })
            .take_while(move |&val| val < end)
    }
}

impl<'a> Iterator for RowBitIter<'a> {
    type Item = u32;

    #[inline(always)]
    fn next(&mut self) -> Option<u32> {
        let mut ret = 0;
        let mut cycles = 0;

        loop {
            let n = self.current.trailing_zeros() + 2;
            self.current >>= n;
            ret += n / 2;

            if 0 == self.current {
                // We counted the sentinel as a cell; fix that
                ret -= 1;

                self.current = prepare_rbi_chunk(unsafe {
                    *self.row.get_unchecked(self.next_chunk)
                });
                self.next_chunk = (self.next_chunk + 1) &
                    (self.row.len() - 1);
                cycles += 1;

                // If the whole row is empty, we can't ever return anything.
                if cycles > self.row.len() {
                    return None;
                }
            } else {
                return Some(ret);
            }
        }
    }
}

/// Wrapper around `RowBitIter` which returns bits 4 at a time.
///
/// Each iterator item is an offset from the *start* of the previous item to
/// the start of the new item, and a bitmask of up to 4 cells from that
/// position, at offsets 0, 2, 1, and 3, in that order.
///
/// Returned bitsets may be incomplete; i.e., a 0 bit does not indicate that
/// the cell is definitely not there; the next iteration may return an item
/// indicating that it does exist.
///
/// Returned bitsets may have arbitrary bits above bit 3 set; these have no
/// particular meaning.
///
/// The iterator will occasionally return completely zero bitsets as well.
///
/// Like `RowBitIter`, this iterator is infinite.
#[derive(Clone, Copy, Debug)]
struct RowNybbleIter<'a> {
    inner: RowBitIter<'a>,
    /// The base offset to apply to the next item. This is used to account for
    /// the actual width of the prior bitset.
    next_off: u32,
}

impl<'a> RowNybbleIter<'a> {
    fn wrap(inner: RowBitIter<'a>) -> Self {
        RowNybbleIter { inner, next_off: 0 }
    }
}

/// Adapt a 4-at-a-time iterator to return column coordinates and bitmasks.
///
/// Zero items are filtered out.
///
/// `end` is only approximate; the iterator may iterate slightly beyond it.
fn scan_through4<IT : Iterator<Item = (u32,u32)>>
    (it: IT, base: i16, end: i16)
     -> impl Iterator<Item = (i32x4,u32)>
{
    it.scan(base, |accum, (step, mask)| {
        *accum += step as i16;
        Some((*accum, mask))
    })
        .take_while(move |&(ix, _)| ix < end)
        .filter(|&(_, mask)| 0 != mask & 0xF)
        .map(|(base, mask)| (i32x4::splat(base as i32) +
                             i32x4::new(0, 2, 1, 3), mask))
}


impl<'a> Iterator for RowNybbleIter<'a> {
    type Item = (u32, u32);

    #[inline(always)]
    fn next(&mut self) -> Option<(u32,u32)> {
        let n = self.inner.current.trailing_zeros();
        let head = self.inner.current >> n;
        // Squash from
        //
        // 76543210
        // -3-2-1-0
        //
        // to
        //
        // 3210
        // 3120
        //
        // i.e, shift 4 to 1 and 6 to 3
        let mut mask = (head | (head >> 3)) as u32;
        let cnt = self.next_off + n / 2;

        self.inner.current >>= n + 8;
        self.next_off = 4;

        // See if this finishes the chunk
        if unlikely!(0 == self.inner.current) {
            // We've consumed the sentinel, remove it
            let lz = head.leading_zeros();
            let mm = 1 << (63 - lz);
            mask ^= mm | (mm >> 3);

            // Fix the stride of the next item
            self.next_off = (63 - lz)/2;

            self.inner.current = prepare_rbi_chunk(unsafe {
                *self.inner.row.get_unchecked(self.inner.next_chunk)
            });
            self.inner.next_chunk = (self.inner.next_chunk + 1) &
                (self.inner.row.len() - 1);
        }

        Some((cnt, mask))
    }
}

/// Similar to `RowNybbleIter`, but only works for composites with a `pitch` of
/// 0 (i.e., single-chunk-wide rows).
///
/// It may return false positives beyond the end of the iterator, and does not
/// detect the end itself.
#[derive(Clone, Copy, Debug)]
struct SingleChunkNybbleIter {
    current: u64,
    next_advance: u32,
}

impl SingleChunkNybbleIter {
    /// Create an iterator starting from the `start`th bit in `chunk`.
    ///
    /// `start` is a bit index as used elsewhere to address chunk bits.
    fn in_row(chunk: Chunk, start: u16) -> Self {
        Chunk::assert_valid_index(start);

        let v = chunk.0 as u64 & 0x555555555555;
        // Rotate the lower `2*CHUNK_WIDTH` bits so that bit 0 is the one that
        // was originally at `start`.
        let v = (v >> start) | (v << CHUNK_WIDTH*2 - start as u32);

        SingleChunkNybbleIter { current: v, next_advance: 0 }
    }
}

impl Iterator for SingleChunkNybbleIter {
    type Item = (u32, u32);

    #[inline(always)]
    fn next(&mut self) -> Option<(u32, u32)> {
        let n = self.current.trailing_zeros();
        let head = self.current >> (n & 63);
        let mask = (head | (head >> 3)) as u32;
        self.current = head >> 8;
        let advance = n/2 + self.next_advance;
        self.next_advance = 4;

        Some((advance, mask))
    }
}

pub const COLLISION_SET_SIZE: usize = 16;
pub type CollisionSet = ArrayVec<[[(i16,i16);2];COLLISION_SET_SIZE]>;

/// The common prefix data shared by all composite objects.
///
/// The first word in a composite is always a `CompositeHeader`.
///
/// This is followed by an array of `Chunk`s of length `1 << rows << pitch`.
///
/// `Chunk`s are grouped into rows, where each row is `1 << pitch` `Chunk`s
/// wide. Given the `A` coordinate of a cell, the cell is found in row
/// `A & ((1 << rows) - 1)`, which allows addressing the rows without prior
/// bounds checks or needing to take `row_offset` into account.
///
/// Each row encodes a neighbourhood bitset for `CHUNK_WIDTH << pitch` cell
/// slots. Given the `B` coordinate of a cell, the cell's neighbourhood is
/// found in slot `B % (CHUNK_WIDTH << pitch)`, i.e., chunk
/// `(B / CHUNK_WIDTH) & ((1 << pitch) - 1)`, chunk slot `B % CHUNK_WIDTH`. As
/// with row access, this means no bounds checking is needed prior to accessing
/// the cell, which is particularly important as the bounds are not known until
/// the chunk is loaded.
///
/// Each row stores the neighbourhood of cells starting at `col_base` and going
/// up to `col_base + (CHUNK_WIDTH << pitch)`. `col_base` varies per column,
/// and is stored as a 14-bit signed integer in every `Chunk`.
///
/// The other 50 bits of a `Chunk` is a bitset of present cells, including
/// their neighbourhood. Given the index of a cell in the chunk, the bit at
/// `index+0` indicates whether that cell is present; `index+1` indicates
/// whether the cell offset by `<1,0>` is present; `index+2`, `<0,1>`, and
/// `index+3`, `<1,1>`. Adjacent cells in the same row have overlapping
/// neighbourhoods, such that each cell only takes two bits of space, except
/// for the final cell in each `Chunk`, which has the left-most two bits of the
/// next cell on the row, though that cell is properly found in the next
/// `Chunk`.
///
/// Every row is sized so that its bitset represents not only the cells on that
/// row, but also the cells on the next row.
///
/// Every row is padded with an empty column at the start and the composite as
/// a whole with an empty row at the top. This means that inspecting the
/// neighbourhood of an out-of-bounds-by-minus-one cell will still reveal
/// populated in-bounds cells. Due to modular indexing, it also ensures that
/// inspecting the last in-bounds cell will have empty bits for its
/// out-of-bounds neighbourhood.
///
/// There is no need to concern with padding after the `Chunk` array, as any
/// non-empty composite will have at least two rows due to the insertion of an
/// empty row at the top.
#[derive(Clone, Copy)]
pub struct CompositeObject<T : Borrow<[i32x4]>>(T);

impl<T : Borrow<[i32x4]>> CompositeObject<T> {
    /// Wraps the given inner array into a `CompositeObject`.
    ///
    /// ## Unsafety
    ///
    /// Bounds checking is only performed in debug builds.
    #[inline(always)]
    pub unsafe fn wrap(inner: T) -> Self {
        {
            let data = inner.borrow();
            debug_assert!(data.len() >= 1);
            let header = CompositeHeader(*data.get_unchecked(0));
            debug_assert!(data.len() >= 1 +
                          (1 << (header.rows() + header.pitch() - 1)));
        }
        CompositeObject(inner)
    }

    #[inline(always)]
    pub fn as_inner(&self) -> &T {
        &self.0
    }

    /// Return a mutable reference to the inner data.
    ///
    /// ## Unsafety
    ///
    /// Mutating the header or reducing the length of the inner value can break
    /// invariants `CompositeObject` relies on.
    #[inline(always)]
    pub unsafe fn as_inner_mut(&mut self) -> &mut T {
        &mut self.0
    }

    #[inline(always)]
    pub fn into_inner(self) -> T {
        self.0
    }

    /// Returns the header of this composite.
    #[inline(always)]
    pub fn header(&self) -> CompositeHeader {
        CompositeHeader(*unsafe { self.0.borrow().get_unchecked(0) })
    }

    /// Returns the slice of data beyond the end of the common composite data.
    #[inline(always)]
    pub fn tail(&self) -> &[i32x4] {
        let s = self.0.borrow();
        let off = 1 + (1 << self.header().rows() << self.header().pitch());
        unsafe {
            slice_unchecked(s, off, s.len())
        }
    }

    #[inline(always)]
    fn chunks(&self) -> &[Chunk] {
        let s = self.0.borrow();
        unsafe {
            slice::from_raw_parts(
                s.as_ptr().offset(1) as *const Chunk,
                (s.len() - 1) * 2)
        }
    }

    /// Returns the physical row index of the given row within the chunk array.
    #[inline(always)]
    fn row_phys_index(&self, row: i16) -> usize {
        (row & ((1 << self.header().rows()) - 1)) as usize
    }

    /// Returns the index of the first chunk belonging to `row` in the chunk
    /// array.
    #[inline(always)]
    fn row_chunk_index(&self, row: i16) -> usize {
        self.row_phys_index(row) << self.header().pitch()
    }

    /// Like `row_chunk_index`, but processes 4 rows at once.
    #[inline(always)]
    fn row_chunk_index4(&self, row: i32x4) -> u32x4 {
        let row = row & i32x4::splat((1 << self.header().rows()) - 1);
        let row: i32x4 = row << (self.header().pitch() as u32);
        u32x4::from_cast(row)
    }

    /// Returns the chunks in the given row.
    #[inline(always)]
    fn row(&self, row: i16) -> &[Chunk] {
        let row_base = self.row_chunk_index(row);
        unsafe {
            self.chunks().get_unchecked(
                row_base..row_base + (1 << self.header().pitch()))
        }
    }

    /// Returns the chunk offset (within the row) and bit (as a multiple of
    /// two) of the given column.
    #[inline(always)]
    fn col_index(&self, col: i16) -> ChunkRowCoord {
        // Since the % operator doesn't work in a sane way for negative
        // integers, convert to unsigned and work with that.
        let col = (col as u16).wrapping_add(32768);
        // NB The optimiser fuses both of these into a single integer multiply
        // and then a couple adds and shifts.
        let d24 = col / 24;
        let bit = (col % 24) << 1;
        let chunk = d24 & ((1 << self.header().pitch()) - 1);
        ChunkRowCoord { chunk, bit }
    }

    /// Like `col_index()`, but translates four columns at a time.
    ///
    /// The chunk index is stored in the high half of each lane, and the bit
    /// index in the low half.
    #[inline(always)]
    fn col_index4(&self, col: i32x4) -> u32x4 {
        // This is the code LLVM generates for `col_index`, translated to SIMD.
        // (u32x4 doesn't have division or modulo operators since SSE has no
        // integer division support.) Variable names are just the registers
        // LLVM happened to choose; don't read anything special into them.
        let esi = u32x4::from_cast(col + i32x4::splat(32768));
        let edx = esi * u32x4::splat(0xaaab);   // imul edx,eax,0xaaab
        let edx: u32x4 = edx >> 0x14;           // shr edx,0x14
        let eax: u32x4 = edx << 3;              // lea eax,[rdx*8+0x0]
        let eax: u32x4 = eax + (eax << 1);      // lea eax,[rax+rax*2]
        let esi = esi - eax;                    // sub esi,eax
        let cl = self.header().pitch();         // Multiple instructions
        let eax = u32x4::splat(1);              // mov eax,0x1
        let eax: u32x4 = eax << (cl as u32);    // shl eax,cl
        // LLVM's code is sightly different for some reason
        let eax = eax - u32x4::splat(1);        // add eax,0xfff
        let eax = eax & edx;                    // and eax,edx

        // esi now holds the bit index, eax the chunk index
        (eax << 16) | (esi << 1)
    }

    /// Load the chunk containing `(row, col)` and return that chunk and the
    /// index of that bit within it.
    #[inline(always)]
    pub fn chunk(&self, row: i16, col: i16) -> (Chunk, u16) {
        let row_base = self.row_chunk_index(row);
        let col_index = self.col_index(col);
        (unsafe {
            *self.chunks().get_unchecked(row_base + col_index.chunk as usize)
        }, col_index.bit)
    }

    /// Return whether the cell at `(a,b)` is populated.
    ///
    /// If `(a,b)` is not in range, the result is unspecified but safe.
    #[inline(always)]
    pub fn is_populated(&self, a: i16, b: i16) -> bool {
        let (chunk, bit) = self.chunk(a, b);
        chunk.is_populated(bit)
    }

    /// Return whether `a` is an in-bounds cell row for this composite.
    ///
    /// The valid range includes the blank row at the beginning but excludes
    /// padding rows at the end.
    #[inline(always)]
    pub fn is_in_a_bound(&self, a: i16) -> bool {
        let header = self.header();
        let row = (a - header.row_offset()) as u32;
        row <= header.row_count() as u32
    }

    /// Return whether `b` is an in-bounds cell column for this composite
    /// within the given chunk.
    #[inline(always)]
    pub fn is_in_b_bound(&self, chunk: Chunk, b: i16) -> bool {
        let off = b as i32 - chunk.col_base() as i32;
        (off as u32) < CHUNK_WIDTH << self.header().pitch()
    }

    /// Returns an iterator over the logical row indices in this composite,
    /// excluding the zeroth row which is always empty.
    pub fn rows(&self) -> impl Iterator<Item = i16> {
        let row_offset = self.header().row_offset() as i32;
        (1..self.header().row_count() as i32 + 1).map(
            move |row| (row + row_offset) as i16)
    }

    /// Returns an iterator over the logical indices of populated columns in
    /// the given logical row.
    pub fn cells_in_row<'a>(&'a self, row: i16)
                            -> impl 'a + Iterator<Item = i16> {
        let (bits, base) = RowBitIter::start_of_row(self, row);
        let end = (base as i16) + (CHUNK_WIDTH << self.header().pitch()) as i16;

        bits.scan_through(base as i16, end)
    }

    fn cells_in_row_rbi(&self, row: i16, first_col: i16, last_col: i16)
                        -> (RowBitIter, i16, i16) {
        let (mut bits, base) = RowBitIter::in_row(self, row);
        let actual_start = max(base as i32, first_col as i32);
        let actual_end = min(
            (base as i32) + ((CHUNK_WIDTH as i32) << self.header().pitch()),
            last_col as i32);

        bits.reset(self.col_index(actual_start as i16));
        (bits, actual_start as i16, actual_end as i16)
    }

    /// Returns an iterator over the logical indices of populated columns in
    /// the given logical row, where the column indices are between `first_col`
    /// (inclusive) and `last_col` (exclusive).
    pub fn cells_in_subrow<'a>(&'a self, row: i16,
                               first_col: i16, last_col: i16)
                               -> impl 'a + Iterator<Item = i16> {
        let (bits, start, end) = self.cells_in_row_rbi(
            row, first_col, last_col);
        bits.scan_through(start, end)
    }

    /// Like `cells_in_subrow`, but return an iterator which produces columns 4
    /// at a time.
    ///
    /// Each quadruplet is accompanied by a mask of which columns actually have
    /// cells.
    ///
    /// This is "approximate" in that it will generally overshoot `last_col`,
    /// potentially including cells beyond it or even out of bounds. In
    /// addition to the iterator, it returns the coordinate of the actual last
    /// column (exclusive) within range.
    pub fn cells4_in_subrow_approx<'a>
        (&'a self, row: i16, first_col: i16, last_col: i16)
         -> (impl 'a + Iterator<Item = (i32x4, u32)>, i16)
    {
        let (bits, start, end) = self.cells_in_row_rbi(
            row, first_col, last_col);
        (scan_through4(RowNybbleIter::wrap(bits), start, end), end)
    }

    /// Like `cells4_in_subrow_approx`, but panics if the pitch of this
    /// composite is not 0.
    pub fn sc_cells4_in_subrow_approx
        (&self, row: i16, first_col: i16, last_col: i16)
         -> (impl Iterator<Item = (i32x4, u32)>, i16)
    {
        debug_assert!(0 == self.header().pitch());
        let row_ix = self.row_chunk_index(row);
        let chunk = unsafe {
            *self.chunks().get_unchecked(row_ix)
        };
        let actual_start = max(chunk.col_base(), first_col);
        let actual_end = min(chunk.col_base() + CHUNK_WIDTH as i16, last_col);
        let bit = self.col_index(actual_start).bit;

        (scan_through4(SingleChunkNybbleIter::in_row(chunk, bit),
                       actual_start, actual_end), actual_end)
    }

    /// Returns an iterator over the logical indices of cells which are
    /// populated in this composite.
    pub fn cells<'a>(&'a self) -> impl 'a + Iterator<Item = (i16,i16)> {
        self.rows().flat_map(
            move |row| self.cells_in_row(row).map(move |col| (row, col)))
    }

    /// Test for a collision with a point particle.
    ///
    /// `relative_pos` is the relative position of the point particle relative
    /// to this object's nominal coordinates, including adjusting for rotation.
    ///
    /// If there is a collision, returns the index of the cell affected.
    #[inline(always)]
    pub fn test_point_collision(&self, relative_pos: Vhs) -> Option<(i16,i16)> {
        // TODO Can maybe be faster via first approximation

        let grid_relative_pos = relative_pos - self.header().offset();
        let (ia, ib) = grid_relative_pos.to_index();

        if !self.is_in_a_bound(ia as i16) { return None; }

        let (chunk, bit) = self.chunk(ia as i16, ib as i16);
        if !self.is_in_b_bound(chunk, ib as i16) ||
            !chunk.is_populated(bit)
        { return None; }

        Some((ia as i16, ib as i16))
    }

    /// Computes an upper bound on the collision span of this object (as per
    /// `CommonObject::rounded_span`, but not yet rounded) from the centre of
    /// gravity.
    pub fn calc_span(&self) -> u32 {
        // Sacrifice a bit of precision to avoid overflow
        // Dropping 4 bits means we can support a max distance of 2**(15+4) =
        // 2**19 (8 screens, far more than we'll ever need) rather than just
        // 2**15 (half a screen). We just need to adjust the final calculation
        // to be a bit more conservative.
        const SHIFT: u32 = 4;

        let header = self.header();
        let offset = header.offset() >> SHIFT;
        let mut max_dist = 0;

        // Search for the cell which is the farthest from the centre of
        // gravity. It could be in any row, so we must scan all the rows.
        for row in self.rows() {
            let mut col_iter = self.cells_in_row(row);
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
        CommonObject::span_of_l2(
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
    ///
    /// `dst` must be an empty `CollisionSet`. It will be populated by this
    /// call. It is not returned as a normal return value as this would incur a
    /// 130-byte memcpy on every call. (In the case where the two objects do
    /// not have any candidate cell pairs, this makes a 4x performance
    /// difference, and a ~2x performance difference when there is only one
    /// candidate cell pair.) This does slow down cases where there are
    /// actually collisions by ~50% worst case, but this is OK since collisions
    /// are rare.
    #[inline(never)]
    pub fn test_composite_collision<R : Borrow<[i32x4]>>(
        &self, dst: &mut CollisionSet, self_obj: CommonObject,
        that: &CompositeObject<R>, that_obj: CommonObject,
        self_inverse_rot_xform: Affine2dH
    ) {
        dst.clear();

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
        let that_span = (that_obj.rounded_span() as i32)
            << ROUNDED_SPAN_SHIFT;
        let that_bounds_self_grid = that_cog_self_grid.repr() +
            i32x4::new(-that_span, -that_span, that_span, that_span);
        let that_bounds_self_grid =
            (that_bounds_self_grid >> CELL_HEX_SHIFT) +
            i32x4::new(-1, -1, 1, 1);

        // Now determine the position of our (0,0) cell within `that`'s grid.
        let that_inverse_rot_xform = Affine2dH::rotate_hex(-that_obj.theta());
        let self_rot_xform = self_inverse_rot_xform.inv_rotate_hex();
        let origin_pos = self_obj.pos().dual() +
            self_rot_xform * self.header().offset().dual();
        let self_origin_that_grid =
            (that_inverse_rot_xform * (origin_pos - that_obj.pos().dual()))
            .single() - that.header().offset();

        // Determine the displacement within `that`'s grid for moving (+1,0)
        // and (0,+1) in our own cell grid.
        //
        // Nominally, this would be
        // ```
        // let grid_displacement = that_inverse_rot_xform * self_rot_xform *
        //    Vhl(CELL_HEX_SIZE, 0, 0, CELL_HEX_SIZE);
        // ```
        //
        // The transform itself stands.
        let grid_displacement_xform = that_inverse_rot_xform * self_rot_xform;
        //
        // However, we don't want to premultiply like this due to precision
        // loss; we'll keep the `AFFINE_POINT` factor in and shift it out
        // later. The transform thus has `AFFINE_POINT+CELL_HEX_SHIFT`
        // significant bits; if we refine these constants later, reconsider
        // whether there's any risk of overflow here. Right now, this leaves
        // 31-20=11 significant bits for cell indices, which is far more than
        // enough.
        debug_assert!(20 == AFFINE_POINT + CELL_HEX_SHIFT as u32 -
                      SUBAFFINE_SHIFT);
        // Let C stand for `CELL_HEX_SIZE`. The Vhl multiply described above
        // computes:
        //
        // ```
        // | a b | | C |   | Ca |           | a b |   | 0 |   | Cb |
        // | c d | | 0 | = | Cc |    and    | c d | = | C | = | Cd |
        // ```
        //
        // In other words, `Vhl(Ca, Cc, Cb, Cd)`. Since `C` is a power of two,
        // we can do the multiply with a simple shift.
        let grid_displacement = Vhl::from_repr(
            grid_displacement_xform.repr().shuf(0, 2, 1, 3) <<
                (CELL_HEX_SHIFT as u32 - SUBAFFINE_SHIFT));

        // Determine the first and last rows to scan, and start tracking the
        // base coordinate for that row.
        let first_row = max(that_bounds_self_grid.extract(0),
                            // First row is always empty
                            self.header().row_offset() as i32 + 1);
        let last_row = min(that_bounds_self_grid.extract(2),
                           self.header().row_offset() as i32 +
                           self.header().row_count() as i32);
        let mut row_zero_off: Vhs =
            grid_displacement.fst() * Vhs(first_row, first_row);

        // Scan the columns of each row for overlapping cells
        // Special-case for pitch=0 composites since they are quite common
        if 0 == self.header().pitch() {
            for row in first_row..last_row + 1 {
                self.test_composite_collision_row_sc(
                    dst, that, grid_displacement,
                    (row_zero_off >> AFFINE_POINT - SUBAFFINE_SHIFT) +
                        self_origin_that_grid,
                    row as i16,
                    that_bounds_self_grid.extract(1) as i16,
                    that_bounds_self_grid.extract(3) as i16);
                row_zero_off = row_zero_off + grid_displacement.fst();
            }
        } else {
            for row in first_row..last_row + 1 {
                self.test_composite_collision_row(
                    dst, that, grid_displacement,
                    (row_zero_off >> AFFINE_POINT - SUBAFFINE_SHIFT) +
                        self_origin_that_grid,
                    row as i16,
                    that_bounds_self_grid.extract(1) as i16,
                    that_bounds_self_grid.extract(3) as i16);
                row_zero_off = row_zero_off + grid_displacement.fst();
            }
        }
    }

    #[inline(always)]
    fn test_composite_collision_row<R : Borrow<[i32x4]>>(
        &self, dst: &mut CollisionSet,
        that: &CompositeObject<R>,
        grid_displacement: Vhl,
        row_zero: Vhs, row: i16,
        first_col: i16, last_col: i16
    ) {
        let (iter, end) = self.cells4_in_subrow_approx(
            row, first_col, last_col);
        for (cols, mask) in iter {
            self.test_composite_collision_col4(
                dst, that, grid_displacement, row_zero, row, cols,
                end, mask);
        }
    }

    #[inline(always)]
    fn test_composite_collision_row_sc<R : Borrow<[i32x4]>>(
        &self, dst: &mut CollisionSet,
        that: &CompositeObject<R>,
        grid_displacement: Vhl,
        row_zero: Vhs, row: i16,
        first_col: i16, last_col: i16
    ) {
        let (iter, end) = self.sc_cells4_in_subrow_approx(
            row, first_col, last_col);
        for (cols, mask) in iter {
            self.test_composite_collision_col4(
                dst, that, grid_displacement, row_zero, row, cols,
                end, mask);
        }
    }

    #[inline(always)]
    fn test_composite_collision_col4<R : Borrow<[i32x4]>>(
        &self, dst: &mut CollisionSet,
        that: &CompositeObject<R>,
        grid_displacement: Vhl,
        row_zero: Vhs,
        row: i16, col: i32x4, last_col: i16,
        mask: u32
    ) {
        let nominal_a = i32x4::splat(row_zero.a()) +
            (i32x4::splat(grid_displacement.snd().a()) * col >>
             AFFINE_POINT - SUBAFFINE_SHIFT);
        let nominal_b = i32x4::splat(row_zero.b()) +
            (i32x4::splat(grid_displacement.snd().b()) * col >>
             AFFINE_POINT - SUBAFFINE_SHIFT);
        let approx_a: i32x4 = nominal_a >> (CELL_HEX_SHIFT as u32);
        let approx_b: i32x4 = nominal_b >> (CELL_HEX_SHIFT as u32);

        let a_in_bounds = approx_a.nsw_between(
            that.header().row_offset() as i32,
            that.header().row_offset() as i32 +
                1 + that.header().row_count() as i32);

        let row_indices = that.row_chunk_index4(approx_a);
        let col_coords = that.col_index4(approx_b);
        let chunk_indices = row_indices + (col_coords >> 16);
        let (chunk0, chunk1, chunk2, chunk3) = unsafe {(
            // AVX512 could do this in one instruction with gather, but any
            // system with AVX512 is so far above the needed spec that there's
            // no reason to bother with it.
            *that.chunks().get_unchecked(chunk_indices.extract(0) as usize),
            *that.chunks().get_unchecked(chunk_indices.extract(1) as usize),
            *that.chunks().get_unchecked(chunk_indices.extract(2) as usize),
            *that.chunks().get_unchecked(chunk_indices.extract(3) as usize)
        )};

        let to_check = mask & a_in_bounds.movemask();
        macro_rules! check_lane {
            ($lane:expr, $chunk:expr) => {{
                if 0 != to_check & (1 << $lane) &&
                    // The column iterator may go out of bounds, so skip if
                    // that happened.
                    col.extract($lane) < last_col as i32
                {
                    let neighbourhood = $chunk.neighbourhood(
                        col_coords.extract($lane) as u16);
                    let b = approx_b.extract($lane);
                    if neighbourhood.any() &&
                        that.is_in_b_bound($chunk, b as i16) &&
                        b as i16 as i32 == b
                    {
                        let hits = neighbourhood.hits(
                            Vhs(nominal_a.extract($lane),
                                nominal_b.extract($lane)));

                        check!(hits, $lane, 0, 0, c00);
                        check!(hits, $lane, 0, 1, c01);
                        check!(hits, $lane, 1, 0, c10);
                        check!(hits, $lane, 1, 1, c11);
                    }
                }
            }}
        }
        macro_rules! check {
            ($hits:expr, $lane:expr, $ao:expr, $bo:expr, $meth:ident) => {
                // No need for bounds checking due to the row and
                // column padding.
                if $hits.$meth() {
                    if dst.try_push([(row as i16, col.extract($lane) as i16),
                                     (approx_a.extract($lane) as i16 + $ao,
                                      approx_b.extract($lane) as i16 + $bo)])
                        .is_err()
                    {
                        return;
                    }
                }
            }
        }
        check_lane!(0, chunk0);
        check_lane!(1, chunk1);
        check_lane!(2, chunk2);
        check_lane!(3, chunk3);
    }
}

impl<T : BorrowMut<[i32x4]>> CompositeObject<T> {
    #[inline(always)]
    fn chunks_mut(&mut self) -> &mut [Chunk] {
        let s = self.0.borrow_mut();
        unsafe {
            slice::from_raw_parts_mut(
                s.as_mut_ptr().offset(1) as *mut Chunk,
                (s.len() - 1) * 2)
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
        // Determine the logical index of the zeroth row. Note that we insert a
        // blank row above.
        let min_row = cells().next().expect(
            "Attempted to make composite with no cells").0 - 1;

        // Make a first pass through the cells to determine the beginning and
        // end logical indices of each row.
        let mut row_spans = SmallVec::<[(i16,i16);64]>::new();
        let mut prev = (i16::MIN, i16::MIN);
        for (a, b) in cells() {
            debug_assert!((a, b) > prev);
            prev = (a, b);

            // Since min_row is offset by 1, this also causes the blank top row
            // to be added implicitly
            let row_ix = (a - min_row) as usize;
            // Add any empty rows as well as this row if needed.
            while row_spans.len() <= row_ix {
                // Minus one since we need an empty column at the start
                row_spans.push((b-1, b-1));
            }

            // Update the maximum column for this row
            row_spans[row_ix].1 = b;
        }

        let num_rows = row_spans.len();

        // Add another empty row at the bottom that's the same as the current
        // last row. This one doesn't get included in the output, but just
        // allows doing things with `windows()`.
        let last_row_span = row_spans[num_rows-1];
        row_spans.push(last_row_span);

        let max_row_span = row_spans.windows(2)
            // Each row needs space for both itself and the cells on the next
            // row.
            .map(|s| 1 + max(s[0].1, s[1].1) - min(s[0].0, s[1].0))
            .max()
            .expect("Empty row_spans");

        fn log2_up(v: i16) -> u8 {
            (16 - (v - 1).leading_zeros()) as u8
        }

        // Set up an empty composite
        let rrows = (num_rows as u16).log2_up();
        base.set_row_offset(min_row as i16);
        base.set_rows(rrows);
        debug_assert!(num_rows <= 65536);
        base.set_row_count((num_rows - 1) as u16);
        let pitch = ((max_row_span as u32 + CHUNK_WIDTH - 1) /
                     CHUNK_WIDTH).log2_up();
        base.set_pitch(pitch);

        let mut dst = CompositeObject({
            let mut data = SmallVec::new();
            data.push(base.0);
            for _ in 0..(1 << rrows-1 << pitch) {
                data.push(i32x4::splat(0));
            }
            data
        });

        // Set all the col_bases for populated rows. We can leave the bases for
        // completely blank rows as zero.
        for (row_ix, spans) in row_spans.windows(2).enumerate() {
            let logical_row = row_ix as i16 + min_row;
            let row_start = dst.row_chunk_index(logical_row);
            for ix in row_start..row_start + (1 << pitch) {
                dst.chunks_mut()[ix].set_col_base(
                    min(spans[0].0, spans[1].0));
            }
        }

        // Finally, populate the bitset
        for (a, b) in cells() {
            macro_rules! poke {
                ($ao:expr, $bo:expr, $meth:ident) => {{
                    let row_start = dst.row_chunk_index(a - $ao);
                    let col_ix = dst.col_index(b - $bo);
                    dst.chunks_mut()[row_start + col_ix.chunk as usize]
                        .$meth(col_ix.bit + 2*$bo);
                }}
            }
            poke!(0, 0, set_col0);
            poke!(0, 1, set_col0);
            poke!(1, 0, set_col1);
            poke!(1, 1, set_col1);
        }

        dst
    }
}

impl<T : Borrow<[i32x4]>> fmt::Debug for CompositeObject<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let header = self.header();

        let mut s = f.debug_struct("CompositeObject");
        s.field("header", &header);

        for row in self.rows() {
            if self.cells_in_row(row).next().is_some() {
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
        for col in self.0.cells_in_row(logical_row) {
            list.entry(&col);
        }
        list.finish()
    }
}

#[cfg(test)]
mod test {
    use std::collections::{BTreeSet, HashSet};
    use std::num::Wrapping;
    use test::{Bencher, black_box};

    use proptest;
    use proptest::strategy::Strategy;

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

    prop_compose! {
        fn arb_composite(max_size: i16)(size in 1..max_size+1)
            (cells in proptest::collection::btree_set(
                (-size..size, -size..size), 1..64),
             a in -65536..65536, b in -65536..65536,
             theta in proptest::num::i16::ANY.prop_map(Wrapping),
             a_offset in -16384..16384,
             b_offset in -16384..16384)
            -> ArbComposite
        {
            let mut ret = ArbComposite {
                common: UnpackedCommonObject {
                    a, b, theta,
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
            let rad = ret.data.calc_span();
            ret.common.set_rounded_span(
                CommonObject::round_span(rad));
            ret
        }
    }

    #[test]
    fn builds_correct_bitset_0neg8() {
        static CELLS: &'static [(i16,i16)] = &[(0,-8)];
        let object = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || CELLS.iter().map(|&v| v))
        };

        assert!(object.is_in_a_bound(0));
        let (chunk, bit) = object.chunk(0, -8);
        assert!(object.is_in_b_bound(chunk, -8));
        assert!(chunk.is_populated(bit));

        let mut it = object.cells();
        assert_eq!(Some((0, -8)), it.next());
        assert_eq!(None, it.next());
    }

    proptest! {
        #[test]
        fn packed_header_preserves_all_fields(
            a_offset in proptest::num::i32::ANY,
            b_offset in proptest::num::i32::ANY,
            row_offset in proptest::num::i16::ANY,
            rows in proptest::num::u8::ANY,
            row_count in proptest::num::u16::ANY,
            pitch in proptest::num::u8::ANY,
            mass in proptest::num::u16::ANY
        ) {
            let orig = UnpackedCompositeHeader {
                a_offset, b_offset, row_offset,
                rows, row_count, pitch, mass,
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
                assert!(object.is_in_a_bound(a),
                        "Cell ({},{}) out of A bound (row_off = {})",
                        a, b, object.header().row_offset());
                let (chunk, bit) = object.chunk(a, b);
                assert!(chunk.is_populated(bit),
                        "Cell ({},{}) not marked populated", a, b);
                assert!(object.is_in_b_bound(chunk, b),
                        "Cell ({},{}) is out of B bound \
                         (col_base = {}, pitch = {})", a, b,
                        chunk.col_base(), 1 << object.header().pitch());

                // Check other neighbourhoods where this cell should be found.
                macro_rules! check {
                    ($ao:expr, $bo:expr, $meth:ident) => {
                        if object.is_in_a_bound(a - $ao) {
                            let (chunk, bit) = object.chunk(a - $ao, b - $bo);
                            if object.is_in_b_bound(chunk, b - $bo) {
                                assert!(chunk.neighbourhood(bit).$meth(),
                                        "Cell ({},{}) not marked present in \
                                         neighbourhood where it should be at \
                                         <{},{}>", a, b, $ao, $bo);
                            }
                        }
                    }
                }
                check!(0, 1, c01);
                check!(1, 0, c10);
                check!(1, 1, c11);
            }

            // Coordinates not covered by `cells` are not populated or are out
            // of bounds
            for a in -128i16..128i16 {
                for b in -128i16..128i16 {
                    if !cells.contains(&(a, b)) {
                        let (chunk, bit) = object.chunk(a, b);
                        let in_a_bound = object.is_in_a_bound(a);
                        let in_b_bound = object.is_in_b_bound(chunk, b);
                        let pop = chunk.is_populated(bit);
                        assert!(!pop || !in_a_bound || !in_b_bound,
                                "Non-cell ({},{}) marked present; \
                                 populated = {}, \
                                 in_a_bound = {}, in_b_bound = {}",
                                a, b, pop, in_a_bound, in_b_bound);
                    }
                }
            }

            // Iterators should return all the cells and no more
            let iterated_cells: BTreeSet<(i16,i16)> = object.cells().collect();
            let not_iterated = cells.difference(&iterated_cells)
                .collect::<Vec<_>>();
            let not_expected = iterated_cells.difference(cells)
                .collect::<Vec<_>>();
            assert!(not_iterated.is_empty(), "Failed to iterate cells {:?}",
                    not_iterated);
            assert!(not_expected.is_empty(), "Iterated NX cells {:?}",
                    not_expected);
        }

        #[test]
        fn simd_col_index_matches_scalar(
            ref composite in arb_composite(40),
            cols in [(-4096i16..4096i16),
                     (-4096i16..4096i16),
                     (-4096i16..4096i16),
                     (-4096i16..4096i16)]
        ) {
            let composite = &composite.data;

            let res = composite.col_index4(
                i32x4::new(cols[0] as i32,
                           cols[1] as i32,
                           cols[2] as i32,
                           cols[3] as i32));
            for i in 0..4 {
                let expected = composite.col_index(cols[i]);
                let actual = ChunkRowCoord {
                    chunk: (res.extract(i as usize) >> 16) as u16,
                    bit: (res.extract(i as usize) & 0xFFFF) as u16,
                };
                assert_eq!(expected, actual);
            }
        }

        #[test]
        fn nybble_iter_same_as_bit_iterator(
            ref c in arb_composite(40)
        ) {
            for row in c.data.rows() {
                let expected = c.data.cells_in_row(row)
                    .collect::<BTreeSet<_>>();
                let min = expected.iter().next().cloned().unwrap_or(0);
                let max = expected.iter().next_back().cloned().unwrap_or(min);

                let mut expected_sloppy = expected.clone();
                expected_sloppy.extend((1..4).map(|v| v+max));

                let actual = c.data.cells4_in_subrow_approx(row, min, max+1)
                    .0
                    .flat_map(|(cols, mask)| (0..4)
                              .filter(move |&lane| 0 != mask & (1 << lane))
                              .map(move |lane| cols.extract(lane) as i16))
                    .collect::<BTreeSet<_>>();

                let not_iterated = expected.difference(&actual)
                    .collect::<Vec<_>>();
                assert!(not_iterated.is_empty(),
                        "Failed to iterate cells {:?}", not_iterated);
                let not_expected = actual.difference(&expected_sloppy)
                    .collect::<Vec<_>>();
                assert!(not_expected.is_empty(),
                        "Iterated unexpected cells {:?}", not_expected);
            }
        }

        #[test]
        fn single_chunk_iter_same_as_bit_iterator(
            ref c in arb_composite(10)
        ) {
            for row in c.data.rows() {
                let expected = c.data.cells_in_row(row)
                    .collect::<BTreeSet<_>>();
                let min = expected.iter().next().cloned().unwrap_or(0);
                let max = expected.iter().next_back().cloned().unwrap_or(min);

                let mut expected_sloppy = expected.clone();
                expected_sloppy.extend((1..4).map(|v| v+max));

                let actual = c.data.sc_cells4_in_subrow_approx(row, min, max+1)
                    .0
                    .flat_map(|(cols, mask)| (0..4)
                              .filter(move |&lane| 0 != mask & (1 << lane))
                              .map(move |lane| cols.extract(lane) as i16))
                    .collect::<BTreeSet<_>>();

                let not_iterated = expected.difference(&actual)
                    .collect::<Vec<_>>();
                assert!(not_iterated.is_empty(),
                        "Failed to iterate cells {:?}", not_iterated);
                let not_expected = actual.difference(&expected_sloppy)
                    .collect::<Vec<_>>();
                assert!(not_expected.is_empty(),
                        "Iterated unexpected cells {:?}", not_expected);
            }
        }
    }

    proptest! {
        #![proptest_config(proptest::test_runner::Config {
            cases: 65536,
            .. proptest::test_runner::Config::default()
        })]

        #[test]
        fn cc_collisions_roughly_correct(
            ref lhs in arb_composite(24),
            ref rhs in arb_composite(24)
        ) {
            // For an exact solution, hexagons would collide somewhere between
            // 2*CELL_L2_EDGE and 2*CELL_L2_VERTEX. Increase the latter by
            // 12/10 to account for various precision loss. The former we
            // reduce to 66% of what it would normally be due to the way the
            // collision test is approximated.
            const AGGRESSIVE_DIST: u32 = CELL_L2_VERTEX as u32 * 2 * 12/10;
            const CONSERVATIVE_DIST: u32 = CELL_L2_EDGE as u32 * 2 * 66/100;

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

            let mut result = CollisionSet::new();
            lhs.data.test_composite_collision(
                &mut result, lhs.common, &rhs.data, rhs.common,
                lhs_inv_rot);
            for item in &result {
                assert!(aggressive_hits.contains(item),
                        "Detected {:?} as a collision, but they should not \
                         have collided. lhs coord = {:?}, rhs coord = {:?}",
                        item, lhs.cell_coord(item[0].0, item[0].1),
                        rhs.cell_coord(item[1].0, item[1].1));
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

    #[bench]
    fn bench_cc_collision_1x1(b: &mut Bencher) {
        static CELLS: &'static [(i16,i16)] = &[(0,0)];

        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || CELLS.iter().map(|&v| v))
        };
        let obj = UnpackedCommonObject {
            rounded_span:
                CommonObject::round_span(composite.calc_span()),
            .. UnpackedCommonObject::default()
        }.pack();
        let xform = Affine2dH::rotate_hex(Wrapping(0));

        b.iter(|| {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision(
                &mut result,
                black_box(obj), black_box(&composite), black_box(obj),
                black_box(xform))
        });
    }

    #[bench]
    fn bench_cc_collision_4x4(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (0..4).flat_map(|r| (0..4).map(move |c| (r, c))))
        };
        let obj = UnpackedCommonObject {
            rounded_span:
                CommonObject::round_span(composite.calc_span()),
            .. UnpackedCommonObject::default()
        }.pack();
        let xform = Affine2dH::rotate_hex(Wrapping(0));


        b.iter(|| {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision(
                &mut result,
                black_box(obj), black_box(&composite), black_box(obj),
                black_box(xform))
        });
    }

    #[bench]
    fn bench_cc_collision_4x4_nonoverlapping(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };
        let obj_a = UnpackedCommonObject {
            a: -5 * CELL_L2_VERTEX, b: -5 * CELL_L2_VERTEX,
            rounded_span:
                CommonObject::round_span(composite.calc_span()),
            .. UnpackedCommonObject::default()
        }.pack();
        let obj_b = UnpackedCommonObject {
            a: 5 * CELL_L2_VERTEX, b: 5 * CELL_L2_VERTEX,
            rounded_span:
                CommonObject::round_span(composite.calc_span()),
            .. UnpackedCommonObject::default()
        }.pack();
        let xform = Affine2dH::rotate_hex(Wrapping(0));


        b.iter(|| {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision(
                &mut result,
                black_box(obj_a), black_box(&composite), black_box(obj_b),
                black_box(xform))
        });
    }

    // Current timing on main test system based on this and more isolated
    // benchmarks:
    //
    // 13ns per col4 call, called 4 times, 52ns
    // 37ns per row call (24ns self), called 4 times, 96ns
    // 48ns before entering the row loop
    // Subtotal is 196ns, reasonably close to the actual time of 203ns
    //
    // Main focus for now then is to continue improving the iterator.
    #[bench]
    fn bench_cc_collision_4x4_noncolliding(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };
        // Artificially inflated span will cause all cells to be checked even
        // though none overlap.
        let obj_a = UnpackedCommonObject {
            a: -5 * CELL_L2_VERTEX, b: -5 * CELL_L2_VERTEX,
            rounded_span: 255,
            .. UnpackedCommonObject::default()
        }.pack();
        let obj_b = UnpackedCommonObject {
            a: 5 * CELL_L2_VERTEX, b: 5 * CELL_L2_VERTEX,
            rounded_span: 255,
            .. UnpackedCommonObject::default()
        }.pack();
        let xform = Affine2dH::rotate_hex(Wrapping(0));

        b.iter(|| {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision(
                &mut result,
                black_box(obj_a), black_box(&composite), black_box(obj_b),
                black_box(xform))
        });
    }

    #[bench]
    fn bench_cc_collision_row(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };

        b.iter(|| {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision_row(
                &mut result, black_box(&composite),
                black_box(Affine2dH::rotate_hex(Wrapping(0)) *
                          Vhl(CELL_HEX_SIZE, 0, 0, CELL_HEX_SIZE)),
                black_box(Vhs(65536, 65536)),
                black_box(0), black_box(-3), black_box(2))
        })
    }

    #[bench]
    fn bench_cc_collision_col4(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };

        b.iter(|| {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision_col4(
                &mut result, black_box(&composite),
                black_box(Affine2dH::rotate_hex(Wrapping(0)) *
                          Vhl(CELL_HEX_SIZE, 0, 0, CELL_HEX_SIZE)),
                black_box(Vhs(65536, 65536)),
                black_box(0),
                black_box(i32x4::new(0, 1, 2, 3)), black_box(3),
                black_box(0xF));
            // Force the function to be evaluated even when inlined. Awkwardly,
            // we can't afford to actually return the `SmallVec` as it turns
            // into a large `memcpy` that completely dominates the function
            // under test.
            result.len()
        })
    }

    #[bench]
    fn bench_composite_col_index(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };
        b.iter(|| composite.col_index(black_box(42)));
    }

    #[bench]
    fn bench_cells_in_subrow_init(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };

        b.iter(|| black_box(&composite).cells4_in_subrow_approx(
            black_box(0), black_box(-4), black_box(4)));
    }

    #[bench]
    fn bench_cells_in_subrow_step(b: &mut Bencher) {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };

        let it = composite.cells4_in_subrow_approx(0, -4, 4).0;

        // The iterator is actually safely cloneable, but closures don't get
        // OBITS yet, so we have to do this unsafely.
        #[inline]
        fn clone_unsafely<T>(t: &T) -> T {
            unsafe { ::std::ptr::read(t as *const T) }
        }

        b.iter(|| black_box(clone_unsafely(&it)).next());
    }

    #[test]
    fn uiaeo() {
        let composite = unsafe {
            CompositeObject::<SmallVec<[i32x4;32]>>::build(
                UnpackedCompositeHeader::default().pack(),
                || (-2..2).flat_map(|r| (-2..2).map(move |c| (r, c))))
        };
        // Artificially inflated span will cause all cells to be checked even
        // though none overlap.
        let obj_a = UnpackedCommonObject {
            a: -5 * CELL_L2_VERTEX, b: -5 * CELL_L2_VERTEX,
            rounded_span: 255,
            .. UnpackedCommonObject::default()
        }.pack();
        let obj_b = UnpackedCommonObject {
            a: 5 * CELL_L2_VERTEX, b: 5 * CELL_L2_VERTEX,
            rounded_span: 255,
            .. UnpackedCommonObject::default()
        }.pack();
        let xform = Affine2dH::rotate_hex(Wrapping(0));

        for _ in 0..1000000 {
            let mut result = CollisionSet::new();
            black_box(&composite).test_composite_collision(
                &mut result,
                black_box(obj_a), black_box(&composite), black_box(obj_b),
                black_box(xform))
        }
    }
}
