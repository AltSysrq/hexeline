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

use simd::*;
use simdext::*;

use physics::common_object::CommonObject;

#[allow(unused_assignments)]
pub fn hilbert_sort(array: &mut [CommonObject]) {
    array.sort_by_key(|o| {
        let pos = o.p;
        // Adapted from
        // https://en.wikipedia.org/wiki/Hilbert_curve#Applications_and_mapping_algorithms
        // For this to be a real implementation, we need to decide what to do
        // with negative coordinates. Also, it needs real tests and such.

        // Drop the lower 8 bits of precision so we need fewer operations
        let mut pos: i32x4 = pos >> 8;
        let mut n = i32x4::splat(1 << 24);
        let mut n1 = n - i32x4::splat(1);
        let mut d = 0u64;

        macro_rules! xy2d {
            ($s:expr) => {
                let r = pos & n;
                d += ((3 * r.extract(0) as u64) ^ r.extract(1) as u64)
                    << ($s * 2);

                if 0 == r.extract(1) {
                    if 1 == r.extract(0) {
                        pos = n1 - pos;
                    }
                    pos = pos.shuf(1, 0, 0, 0);
                }

                n = n >> 1;
                n1 = n1 >> 1;
            }
        }

        xy2d!(24);
        xy2d!(23);
        xy2d!(22);
        xy2d!(21);
        xy2d!(20);
        xy2d!(19);
        xy2d!(18);
        xy2d!(17);
        xy2d!(16);
        xy2d!(15);
        xy2d!(14);
        xy2d!(13);
        xy2d!(12);
        xy2d!(11);
        xy2d!(10);
        xy2d!( 9);
        xy2d!( 8);
        xy2d!( 7);
        xy2d!( 6);
        xy2d!( 5);
        xy2d!( 4);
        xy2d!( 3);
        xy2d!( 2);
        xy2d!( 1);
        xy2d!( 0);
        d
    });
}

#[cfg(test)]
mod test {
    use test::Bencher;

    use physics::common_object::UnpackedCommonObject;
    use super::*;

    #[bench]
    fn bench_hilbert_sort(b: &mut Bencher) {
        let mut data: Vec<CommonObject> =
            (0..10000).into_iter().map(|i| UnpackedCommonObject {
                a: i * 256, .. UnpackedCommonObject::default()
            }.pack()).collect();
        b.iter(|| hilbert_sort(&mut data));
    }
}
