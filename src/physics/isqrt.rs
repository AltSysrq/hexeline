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

use std::i32;

#[inline]
pub fn isqrt32_up(n: u32) -> u32 {
    // The algorithm below will work with n > i32::MAX, but `result_squared`
    // would need to be a u64, and we need to deal with `shift == 32`.
    debug_assert!(n <= i32::MAX as u32);

    // Adapted from
    // https://en.wikipedia.org/wiki/Integer_square_root#Using_bitwise_operations
    // The algorithm there rounds _down_, but we want to round _up_.
    let mut shift = 0;

    // Base case: n <= 1 -> n
    // Reduce larger values to the base case by shifting
    while n >> shift > 1 { shift += 2; }
    // We can solve the base case directly, then build up from there
    let mut result = n >> shift;
    // Compute and maintain `result*result` without using multiplication
    let mut result_squared = result; // 0 or 1 squared

    // Walk the shift back, adjusting `result` as necessary.
    while shift > 0 {
        shift -= 2;
        result <<= 1;
        result_squared <<= 2;
        // If (result+1)**2 is <= n, increase result by 1. Note that this will
        // compute the _floor_ of the square root, which we fix up after the
        // loop.
        // (result+1)*(result+1) =
        // result*result + 2*result + 1
        if result_squared + 2*result + 1 <= n >> shift {
            result_squared += 2 * result + 1;
            result += 1;
        }
    }

    // We found the greatest integer where `result**2 <= n`, but we want the
    // least integer where `result**2 >= n`, so adjust as needed.
    if result_squared < n {
        result + 1
    } else {
        result
    }
}

#[cfg(test)]
mod test {
    use test::{Bencher, black_box};

    use super::*;

    #[test]
    fn isqrt32_up_20() {
        assert_eq!(5, isqrt32_up(20));
    }

    #[test]
    fn isqrt32_up_536870912() {
        assert_eq!(23171, isqrt32_up(536870912));
    }

    proptest! {
        #[test]
        fn test_isqrt32_up(input in 0u32..(1 << 31)) {
            let expected = (input as f64).sqrt().ceil() as u32;
            let actual = isqrt32_up(input);
            assert_eq!(expected, actual);
        }
    }

    #[bench]
    fn bench_isqrt32_up(b: &mut Bencher) {
        b.iter(|| isqrt32_up(black_box(123456789)));
    }

    #[bench]
    fn bench_isqrt32_via_float(b: &mut Bencher) {
        b.iter(|| (black_box(123456789) as f64).sqrt().ceil() as u32)
    }
}
