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

use std::mem;

use simd::*;

pub trait SimdExt {
    fn max(self, that: Self) -> Self;
    fn min(self, that: Self) -> Self;
    fn movemask(self) -> u32;
}

impl SimdExt for i32x4 {
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
}
