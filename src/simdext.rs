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

use std::num::Wrapping;
use std::mem;

use simd::*;

pub trait SimdExt {
    type Lane;
    type ULane;

    /// Return a SIMD value containing the maximum value of each lane.
    fn max(self, that: Self) -> Self;
    /// Return a SIMD value containing the minimum value of each lane.
    fn min(self, that: Self) -> Self;
    /// Return an n-bit value (starting with bit 0) containing the sign bit of
    /// each lane in this SIMD value.
    fn movemask(self) -> u32;
    /// Compute the absolute value of each lane in this value.
    fn abs(self) -> Self;

    /// Return the sum of the first two lanes.
    fn hsum_2(self) -> Self::Lane;
    /// Return the sum of the first three lanes.
    fn hsum_3(self) -> Self::Lane;

    /// Return the 2D L1 distance between the two values.
    ///
    /// In cartesian space, this describes a diamond. In hexagonal space, it is
    /// a sheared diamond.
    #[allow(non_snake_case)]
    fn dist_2L1(self, that: Self) -> Self::ULane;
    /// Return the 3D L1 distance between the two values.
    ///
    /// In hexagonal space, this describes a diamond.
    #[allow(non_snake_case)]
    fn dist_3L1(self, that: Self) -> Self::ULane;
    /// Return the 2D L2 distance squared between the two values.
    ///
    /// In cartesian space, this describes a circle; in hexagonal space, an
    /// oval.
    #[allow(non_snake_case)]
    fn dist_2L2_squared(self, that: Self) -> Self::ULane;
    /// Return the 3D L2 distance squared between the two values.
    ///
    /// In hexagonal space, this describes a circle.
    #[allow(non_snake_case)]
    fn dist_3L2_squared(self, that: Self) -> Self::ULane;
    /// Return the 2D L-infinity distance between the two values.
    ///
    /// In cartesian space, this describes a square; in hexagonal space, a
    /// rhombus.
    #[allow(non_snake_case)]
    fn dist_2Linf(self, that: Self) -> Self::ULane;
    /// Return the 3D L-infinity distance between the two values.
    ///
    /// In hexagonal space, this describes a regular hexagon rotated 30 degrees
    /// relative to the normal hexagonal grid.
    #[allow(non_snake_case)]
    fn dist_3Linf(self, that: Self) -> Self::ULane;
}

pub trait SimdExt4 {
    /// Shuffle self using the given immediates.
    fn shuf(self, a: u32, b: u32, c: u32, d: u32) -> Self;
}

impl SimdExt for i32x4 {
    type Lane = i32;
    type ULane = u32;

    #[inline(always)]
    fn max(self, that: i32x4) -> i32x4 {
        i32x4_max(self, that)
    }

    #[inline(always)]
    fn min(self, that: i32x4) -> i32x4 {
        i32x4_min(self, that)
    }

    #[inline(always)]
    fn movemask(self) -> u32 {
        i32x4_movemask(self)
    }

    #[inline(always)]
    fn abs(self) -> i32x4 {
        i32x4_abs(self)
    }

    #[inline(always)]
    fn hsum_2(self) -> i32 {
        (Wrapping(self.extract(0)) +
         Wrapping(self.extract(1))).0
    }

    #[inline(always)]
    fn hsum_3(self) -> i32 {
        (Wrapping(self.extract(0)) +
         Wrapping(self.extract(1)) +
         Wrapping(self.extract(2))).0
    }

    #[inline]
    fn dist_2L1(self, that: i32x4) -> u32 {
        (self - that).abs().hsum_2() as u32
    }

    #[inline]
    fn dist_3L1(self, that: i32x4) -> u32 {
        (self - that).abs().hsum_3() as u32
    }

    #[inline]
    fn dist_2L2_squared(self, that: i32x4) -> u32 {
        let diff = self - that;
        (diff * diff).hsum_2() as u32
    }

    #[inline]
    fn dist_3L2_squared(self, that: i32x4) -> u32 {
        let diff = self - that;
        (diff * diff).hsum_3() as u32
    }

    #[inline]
    fn dist_2Linf(self, that: i32x4) -> u32 {
        use std::cmp::max;

        let diff = (self - that).abs();
        max(diff.extract(0) as u32, diff.extract(1) as u32)
    }

    #[inline]
    fn dist_3Linf(self, that: i32x4) -> u32 {
        use std::cmp::max;

        let diff = (self - that).abs();
        max(max(diff.extract(0) as u32, diff.extract(1) as u32),
            diff.extract(2) as u32)
    }
}

impl SimdExt4 for i32x4 {
    #[inline(always)]
    fn shuf(self, a: u32, b: u32, c: u32, d: u32) -> i32x4 {
        i32x4::new(self.extract(a), self.extract(b),
                   self.extract(c), self.extract(d))
    }
}

#[cfg(target_feature = "sse4.1")]
#[inline(always)]
fn i32x4_max(a: i32x4, b: i32x4) -> i32x4 {
    use simd::x86::sse4_1::Sse41I32x4;
    Sse41I32x4::max(a, b)
}

#[cfg(target_feature = "sse4.1")]
#[inline(always)]
fn i32x4_min(a: i32x4, b: i32x4) -> i32x4 {
    use simd::x86::sse4_1::Sse41I32x4;
    Sse41I32x4::min(a, b)
}

#[cfg(not(target_feature = "sse4.1"))]
#[inline(always)]
fn i32x4_max(a: i32x4, b: i32x4) -> i32x4 {
    bool32ix4::from_repr((a - b) >> 31).select(b, a)
}

#[cfg(not(target_feature = "sse4.1"))]
#[inline(always)]
fn i32x4_min(a: i32x4, b: i32x4) -> i32x4 {
    bool32ix4::from_repr((a - b) >> 31).select(a, b)
}

#[cfg(target_feature = "sse")]
#[inline(always)]
fn i32x4_movemask(a: i32x4) -> u32 {
    extern "platform-intrinsic"{
        fn x86_mm_movemask_ps(a: f32x4) -> i32;
    }
    unsafe {
        x86_mm_movemask_ps(mem::transmute(a)) as u32
    }
}

#[cfg(not(target_feature = "sse"))]
#[inline(always)]
fn i32x4_movemask(a: i32x4) -> u32 {
    let bits = bool32ix4::from_repr(a >> 31);
    ((bits.extract(0) as u32) << 0) |
    ((bits.extract(1) as u32) << 1) |
    ((bits.extract(2) as u32) << 2) |
    ((bits.extract(3) as u32) << 3)
}

#[cfg(target_feature = "ssse3")]
fn i32x4_abs(a: i32x4) -> i32x4 {
    use simd::x86::ssse3::Ssse3I32x4;
    Ssse3I32x4::abs(a)
}

#[cfg(not(target_feature = "ssse3"))]
fn i32x4_abs(a: i32x4) -> i32x4 {
    bool32ix4::from_repr(a >> 31).select(
        i32x4::splat(0) - a, a)
}

#[cfg(test)]
mod test {
    use simd::*;
    use super::*;

    #[test]
    fn test_i32x4_max() {
        let a = i32x4::new(0, 1, -1, -2);
        let b = i32x4::new(0, 0, 0, -1);
        let v = a.max(b);

        assert_eq!(0, v.extract(0));
        assert_eq!(1, v.extract(1));
        assert_eq!(0, v.extract(2));
        assert_eq!(-1, v.extract(3));
    }

    #[test]
    fn test_i32x4_min() {
        let a = i32x4::new(0, 1, -1, -2);
        let b = i32x4::new(0, 0, 0, -1);
        let v = a.min(b);

        assert_eq!(0, v.extract(0));
        assert_eq!(0, v.extract(1));
        assert_eq!(-1, v.extract(2));
        assert_eq!(-2, v.extract(3));
    }

    #[test]
    fn test_i32x4_movemask() {
        let a = i32x4::new(0, -1, 0x7FFFFFFF, -0x8000000);
        assert_eq!(10, a.movemask());
    }

    #[test]
    fn test_i32x4_abs() {
        let a = i32x4::new(0, 42, -1, 0x7FFFFFFF);
        let v = a.abs();
        assert_eq!(0, v.extract(0));
        assert_eq!(42, v.extract(1));
        assert_eq!(1, v.extract(2));
        assert_eq!(0x7FFFFFFF, v.extract(3));
    }

    #[test]
    fn test_i32x4_hsum() {
        let a = i32x4::new(1, -2, 3, -4);
        assert_eq!(-1, a.hsum_2());
        assert_eq!(2, a.hsum_3());
    }

    #[test]
    fn test_i32x4_dist() {
        let a = i32x4::new(1, 2, 3, 4);
        let b = i32x4::new(-1, 5, -10, 20);

        assert_eq!(5, a.dist_2L1(b));
        assert_eq!(18, a.dist_3L1(b));
        assert_eq!(13, a.dist_2L2_squared(b));
        assert_eq!(182, a.dist_3L2_squared(b));
        assert_eq!(3, a.dist_2Linf(b));
        assert_eq!(13, a.dist_3Linf(b));
    }
}
