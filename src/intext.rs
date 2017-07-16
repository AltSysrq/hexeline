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

pub trait UIntExt {
    /// Return the least integer `n` such that `n*n >= self`.
    fn sqrt_up(self) -> Self;
}

impl UIntExt for u32 {
    #[inline(always)]
    fn sqrt_up(self) -> Self {
        // Flip through floating-point not because it is easy, but because it
        // is *faster* than implementing it with integers in software.
        //
        // See commit ff2ecb93a4bec6d20938510184a20a34434edb34 for such a
        // software implementation. It was in all cases inferior to this
        // method. Some very rough numbers:
        //
        // - AMD64 SSE4.1: FP sqrt 8x faster
        // - AMD64 SSE2: FP sqrt 3x faster
        // - ARM7h: FP sqrt 2x faster
        //
        // All `u32`s can be accurately represented in an f64, so there is no
        // loss of precision by doing this either.
        (self as f64).sqrt().ceil() as u32
    }
}

#[cfg(test)]
mod test {
    use std::u32;
    use test::{Bencher, black_box};

    use proptest;

    use super::*;

    fn u32_sqrt_up_case(i: u32) {
        let s = i.sqrt_up() as u64;
        assert!(s*s >= i as u64 && (0 == s || (s-1)*(s-1) < i as u64),
                "sqrt({}) => {}", i, s);
    }

    #[test]
    fn u32_sqrt_up_smoke_check() {
        for &i in &[0, 1, 2, 3, 4, 65535, 65536, 65537,
                    (1 << 31) - 1, 1 << 31, (1 << 31) + 1,
                    u32::MAX] {
            u32_sqrt_up_case(i);
        }
    }

    proptest! {
        #[test]
        fn u32_sqrt_up_prop(i in proptest::num::u32::ANY) {
            u32_sqrt_up_case(i);
        }
    }

    #[bench]
    fn bench_u32_sqrt_up(b: &mut Bencher) {
        b.iter(|| black_box(42u32).sqrt_up());
    }
}
