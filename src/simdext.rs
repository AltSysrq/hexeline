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

use simd::*;

pub trait SimdExt {
    type Lane;
    type ULane;

    /// Return a SIMD value containing the maximum value of each lane, assuming
    /// that subtracting corresponding lanes never overflows.
    fn nsw_max(self, that: Self) -> Self;
    /// Return a SIMD value containing the minimum value of each lane, assuming
    /// that subtracting corresponding lanes never overflows.
    fn nsw_min(self, that: Self) -> Self;
    /// Clamp the first 3 lanes of this value to be between lower and upper
    /// (both inclusive), and leave the final lane alone. The final lane must
    /// have the minimum value for the type in `lower` and the maximum value
    /// for the type in `upper`.
    fn nsw_clamp3(self, lower: Self, upper: Self) -> Self;
    /// Return an n-bit value (starting with bit 0) containing the sign bit of
    /// each lane in this SIMD value.
    fn movemask(self) -> u32;
    /// Return whether the sign bit is set in any lane.
    fn any_sign_bit(self) -> bool;
    /// Compute the absolute value of each lane in this value.
    fn abs(self) -> Self;
    /// Return a vector with the sign bit set for each lane with falls within
    /// the given half-open range.
    ///
    /// The other bits in the result have no particular value. It is intended
    /// to either be sign-extended into a mask or passed to `movemask()`.
    ///
    /// Addition and subtraction on the inputs are assumed to not overflow.
    fn nsw_between(self, lo: Self::Lane, hi: Self::Lane) -> Self;

    /// Return the sum of the first two lanes.
    fn hsum_2(self) -> Self::Lane;
    /// Return the sum of the first three lanes.
    fn hsum_3(self) -> Self::Lane;

    /// Compute self * that / 2**point without intermediate overflow or loss of
    /// precision.
    fn mulfp(self, small: Self, point: u32) -> Self;

    /// Compute the double-width product of corresponding even-numbered lanes
    /// in the two vectors. Store the lower half of each product in the
    /// corresponding even-numbered lanes and the upper half of each product in
    /// the following odd-numbered lanes.
    fn muld(self, rhs: Self) -> Self;

    /// Interpret this vector as little-endian double words and perform a
    /// arithmetic right shift by the given amount.
    fn double_shar(self, n: u32) -> Self;

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

    /// Produce a single vector by blending the two vectors together.
    ///
    /// Each index can be 00, 01, 02, or 03 to select from the left vector, or
    /// 10, 11, 12, or 13 to select from the second vector.
    fn blend(self, right: Self, a: u32, b: u32, c: u32, d:u32) -> Self;
}

impl SimdExt for i32x4 {
    type Lane = i32;
    type ULane = u32;

    #[inline(always)]
    fn nsw_max(self, that: i32x4) -> i32x4 {
        i32x4_nsw_max(self, that)
    }

    #[inline(always)]
    fn nsw_min(self, that: i32x4) -> i32x4 {
        i32x4_nsw_min(self, that)
    }

    #[inline(always)]
    fn nsw_clamp3(self, lower: i32x4, upper: i32x4) -> i32x4 {
        i32x4_nsw_clamp3(self, lower, upper)
    }

    #[inline(always)]
    fn movemask(self) -> u32 {
        i32x4_movemask(self)
    }

    #[inline(always)]
    fn any_sign_bit(self) -> bool {
        i32x4_any_sign_bit(self)
    }

    #[inline(always)]
    fn abs(self) -> i32x4 {
        i32x4_abs(self)
    }

    #[inline(always)]
    fn nsw_between(self, lo: i32, hi: i32) -> i32x4 {
        let cmp_lo = i32x4::splat(lo-1) - self;
        let cmp_hi = self - i32x4::splat(hi);
        cmp_lo & cmp_hi
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

    #[inline(always)]
    fn mulfp(self, rhs: i32x4, point: u32) -> i32x4 {
        let p02 = self.muld(rhs);
        let p13 = self.shuf(1, 0, 3, 2).muld(rhs.shuf(1, 0, 3, 2));
        let p02 = p02.double_shar(point);
        let p13 = p13.double_shar(point);
        p02.blend(p13, 00, 10, 02, 12)
    }

    #[inline(always)]
    fn muld(self, rhs: i32x4) -> i32x4 {
        i32x4_muld(self, rhs)
    }

    #[inline(always)]
    fn double_shar(self, n: u32) -> i32x4 {
        i32x4_double_shar(self, n)
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

    #[inline(always)]
    fn blend(self, other: Self, a: u32, b: u32, c: u32, d: u32) -> Self {
        macro_rules! which {
            ($v:expr) => { if $v < 10 {
                self.extract($v)
            } else {
                other.extract($v - 10)
            } }
        }

        i32x4::new(which!(a), which!(b), which!(c), which!(d))
    }
}

#[cfg(target_feature = "sse4.1")]
#[inline(always)]
fn i32x4_nsw_max(a: i32x4, b: i32x4) -> i32x4 {
    use simd::x86::sse4_1::Sse41I32x4;
    Sse41I32x4::max(a, b)
}

#[cfg(target_feature = "sse4.1")]
#[inline(always)]
fn i32x4_nsw_min(a: i32x4, b: i32x4) -> i32x4 {
    use simd::x86::sse4_1::Sse41I32x4;
    Sse41I32x4::min(a, b)
}

#[cfg(all(not(target_feature = "sse4.1"), not(target_feature = "neon")))]
#[inline(always)]
fn i32x4_nsw_max(a: i32x4, b: i32x4) -> i32x4 {
    bool32ix4::from_repr((a - b) >> 31).select(b, a)
}

#[cfg(all(not(target_feature = "sse4.1"), not(target_feature = "neon")))]
#[inline(always)]
fn i32x4_nsw_min(a: i32x4, b: i32x4) -> i32x4 {
    bool32ix4::from_repr((a - b) >> 31).select(a, b)
}

#[cfg(target_feature = "neon")]
#[inline(always)]
fn i32x4_nsw_max(a: i32x4, b: i32x4) -> i32x4 {
    extern "platform-intrinsic" {
        fn arm_vmaxq_s32(a: i32x4, b: i32x4) -> i32x4;
    }
    unsafe {
        arm_vmaxq_s32(a, b)
    }
}

#[cfg(target_feature = "neon")]
#[inline(always)]
fn i32x4_nsw_min(a: i32x4, b: i32x4) -> i32x4 {
    extern "platform-intrinsic" {
        fn arm_vminq_s32(a: i32x4, b: i32x4) -> i32x4;
    }
    unsafe {
        arm_vminq_s32(a, b)
    }
}

#[cfg(any(target_feature = "sse4.1", target_feature = "neon"))]
#[inline(always)]
fn i32x4_nsw_clamp3(a: i32x4, lower: i32x4, upper: i32x4) -> i32x4 {
    a.nsw_max(lower).nsw_min(upper)
}

#[cfg(not(any(target_feature = "sse4.1", target_feature = "neon")))]
#[inline(always)]
fn i32x4_nsw_clamp3(a: i32x4, lower: i32x4, upper: i32x4) -> i32x4 {
    let clamp3 = a.nsw_max(lower).nsw_min(upper);
    i32x4::new(clamp3.extract(0), clamp3.extract(1), clamp3.extract(2),
               a.extract(3))
}

#[cfg(target_feature = "sse")]
#[inline(always)]
fn i32x4_movemask(a: i32x4) -> u32 {
    extern "platform-intrinsic"{
        fn x86_mm_movemask_ps(a: f32x4) -> i32;
    }
    unsafe {
        use std::mem;
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

#[cfg(target_feature = "sse")]
#[inline(always)]
fn i32x4_any_sign_bit(a: i32x4) -> bool {
    0 != a.movemask()
}

#[cfg(not(target_feature = "sse"))]
#[inline(always)]
fn i32x4_any_sign_bit(a: i32x4) -> bool {
    bool32ix4::from_repr(a >> 31).any()
}

#[cfg(target_feature = "ssse3")]
#[inline(always)]
fn i32x4_abs(a: i32x4) -> i32x4 {
    use simd::x86::ssse3::Ssse3I32x4;
    Ssse3I32x4::abs(a)
}

#[cfg(not(target_feature = "ssse3"))]
#[inline(always)]
fn i32x4_abs(a: i32x4) -> i32x4 {
    bool32ix4::from_repr(a >> 31).select(
        i32x4::splat(0) - a, a)
}

#[cfg(target_feature = "sse2")]
#[inline(always)]
fn i32x4_muld(a: i32x4, b: i32x4) -> i32x4 {
    use std::mem::transmute;
    use simd::x86::sse2::i64x2;

    extern "platform-intrinsic"{
        // Weirdly, this definition differs from the Intel spec, though the
        // distinction is largely moot.
        fn x86_mm_mul_epi32(a: i32x4, b: i32x4) -> i64x2;
    }

    unsafe {
        transmute(x86_mm_mul_epi32(a, b))
    }
}

#[cfg(not(target_feature = "sse2"))]
#[inline(always)]
fn i32x4_muld(a: i32x4, b: i32x4) -> i32x4 {
    // Not all platforms have i64x2, but LLVM can emulate it when not.
    // Unfortunately the simd crate doesn't expose this, so define by hand.
    extern "platform-intrinsic" {
        fn simd_mul<T>(a: T, b: T) -> T;
        fn simd_extract<T, U>(a: T, ix: u32) -> U;
    }
    #[repr(simd)]
    struct I64x2(i64, i64);
    #[inline(always)]
    fn extract(v: I64x2, ix: u32) -> i64 {
        simd_extract(v, ix)
    }

    let a = I64x2(a.extract(0) as i64, a.extract(2) as i64);
    let b = I64x2(b.extract(0) as i64, b.extract(2) as i64);
    let r = unsafe { simd_mul(a, b) };
    i32x4::new(extract(r,0) as i32, (extract(r,0) >> 32) as i32,
               extract(r,1) as i32, (extract(r,1) >> 32) as i32)
}

#[cfg(target_feature = "sse2")]
#[inline(always)]
fn i32x4_double_shar(a: i32x4, n: u32) -> i32x4 {
    use std::mem::transmute;
    use simd::x86::sse2::i64x2;

    unsafe {
        let a: i64x2 = transmute(a);
        let a = a >> n;
        transmute(a)
    }
}

#[cfg(not(target_feature = "sse2"))]
#[inline(always)]
fn i32x4_double_shar(a: i32x4, n: u32) -> i32x4 {
    let base: i32x4 = (a.to_u32() >> n).to_i32();
    let upper: i32x4 = a << 32 - n;
    base | i32x4::new(upper.extract(1), 0, upper.extract(3), 0)
}

#[cfg(test)]
mod test {
    use std::i32;
    use test::{Bencher, black_box};

    use simd::*;

    use proptest;

    use super::*;

    #[test]
    fn test_i32x4_max() {
        let a = i32x4::new(0, 1, -1, -2);
        let b = i32x4::new(0, 0, 0, -1);
        let v = a.nsw_max(b);

        assert_eq!(0, v.extract(0));
        assert_eq!(1, v.extract(1));
        assert_eq!(0, v.extract(2));
        assert_eq!(-1, v.extract(3));
    }

    #[test]
    fn test_i32x4_min() {
        let a = i32x4::new(0, 1, -1, -2);
        let b = i32x4::new(0, 0, 0, -1);
        let v = a.nsw_min(b);

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
    fn test_i32x4_nsw_between() {
        assert_eq!(0xF, i32x4::splat(0).nsw_between(-1, 1).movemask());
        assert_eq!(0xF, i32x4::splat(0).nsw_between(0, 1).movemask());
        assert_eq!(0xF, i32x4::splat(-1).nsw_between(-1, 1).movemask());
        assert_eq!(0x0, i32x4::splat(-1).nsw_between(0, 1).movemask());
        assert_eq!(0x0, i32x4::splat(1).nsw_between(-1, 1).movemask());
        assert_eq!(0x0, i32x4::splat(2).nsw_between(-1, 1).movemask());
    }

    #[test]
    fn test_i32x4_hsum() {
        let a = i32x4::new(1, -2, 3, -4);
        assert_eq!(-1, a.hsum_2());
        assert_eq!(2, a.hsum_3());
    }

    #[test]
    fn test_i32x4_mulfp() {
        let mut fibs = vec![0i32, 0i32, 1i32, -1i32];
        loop {
            let a = fibs[fibs.len() - 4] as u64;
            let b = fibs[fibs.len() - 2] as u64;
            let f = a + b;
            if f <= i32::MAX as i64 as u64 {
                fibs.push(f as i32);
                fibs.push(-(f as i32));
            } else {
                break;
            }
        }
        fibs.push(32767);
        fibs.push(32768);

        for &lhs in &fibs {
            for &rhs in &fibs {
                for point in 1..16 {
                    if rhs.abs() > 32768 { continue; }

                    let expected = (lhs as i64) * (rhs as i64) >> point;
                    if expected > i32::MAX as i64 ||
                        expected < i32::MIN as i64 { continue; }

                    let actual = i32x4::splat(lhs).mulfp(
                        i32x4::splat(rhs), point);
                    assert!(actual.extract(0) as i64 == expected,
                            "Expected {} * {} >> {} = {}, got {}",
                            lhs, rhs, point, expected, actual.extract(0));
                    assert_eq!(actual.extract(0), actual.extract(1));
                    assert_eq!(actual.extract(0), actual.extract(2));
                    assert_eq!(actual.extract(0), actual.extract(3));
                }
            }
        }
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

    proptest! {
        #[test]
        fn test_i32x4_muld(a in proptest::num::i32::ANY,
                           b in proptest::num::i32::ANY) {
            let p = (a as u64).wrapping_mul(b as u64);
            let expected_lo = p as i32;
            let expected_hi = (p >> 32) as i32;

            let actual = i32x4::new(a, 99, a, -99).muld(
                i32x4::new(b, 88, b, -88));
            assert_eq!(expected_lo, actual.extract(0));
            assert_eq!(expected_lo, actual.extract(2));
            assert_eq!(expected_hi, actual.extract(1));
            assert_eq!(expected_hi, actual.extract(3));
        }
    }

    macro_rules! bench_dist {
        ($test_name:ident, $dist_f:ident) => {
            #[bench]
            fn $test_name(b: &mut Bencher) {
                b.iter(|| {
                    let a = i32x4::new(1, 2, 3, 4);
                    let b = i32x4::new(5, 6, 7, 8);
                    // 10 times
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b)) +
                    black_box(a).$dist_f(black_box(b))
                })
            }
        }
    }

    bench_dist!(bench_bounding_shear_diamond, dist_2L1);
    bench_dist!(bench_bounding_diamond, dist_3L1);
    bench_dist!(bench_bounding_oval, dist_2L2_squared);
    bench_dist!(bench_bounding_circle, dist_3L2_squared);
    bench_dist!(bench_bounding_rhombus, dist_2Linf);
    bench_dist!(bench_bounding_hexagon, dist_3Linf);
}
