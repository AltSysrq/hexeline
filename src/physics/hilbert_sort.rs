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

#[inline(never)]
#[allow(unused_assignments)]
fn xy_to_hilbert(pos: i32x4) -> u64 {
    /*
    On https://en.wikipedia.org/wiki/Hilbert_curve#Applications_and_mapping_algorithms
    we have an algorithm that's essentially

      let mut d = 0;
      for bit in 24...0 {
        let r = (xy & (1 << bit)) >> bit;
        d += ((3 * r.extract(0)) ^ r.extract(1)) << bit << bit;
        if 0 == r.extract(1) {
          if 0 != r.extract(0) {
            xy = i32x4::splat((1 << bit) - 1) - xy;
          }
          xy = xy.shuf(1, 0);
        }
      }
      d

    The main unfortunate things about this is requiring 32 iterations, as well
    as two branches per iteration. The first thing to do is eliminate the branches.

      let r_mask = r << 31 >> 31;
      let inverted = i32x4::splat((1 << bit) - 1) - xy;
      let alt_inverted = r_mask.shut(0, 0, 0, 0).select(xy, inverted);
      let alt = alt_inverted = alt.shuf(1, 0, 0, 0);
      xy = r_mask.shuf(1, 1, 1, 1).select(alt, xy);

    Next, note that we're currently only utilising 2 of 4 lanes. If we could
    use all 4, we could do two operations concurrently, making the algorithm
    almost twice as fast. But if we can split the problem in two, we can split
    it in eight and use all 16 lanes of an u8x16.

    The first obstacle to doing this is how the results are combined. This is
    actually fairly easy; each bit of input uniquely contributes to two bits of
    output, so there is no carry between the lanes. Combining the outputs is
    thus just a matter of pasting the bits together. If we put 4 input bits
    into each lane, we get 8 bits of output from each, and if we order the
    lanes properly, the output can be obtained for free with a transmute.

    The other big obstacle is that the process for each bit appears to depend
    on all the bits preceeding it. The first thing to note is that there are
    two distinct operations performed on the input coordinates:

    - Swap X and Y
    - Invert all bits below `bit` in both X and Y

    These can *almost* be performed independently. The problem is that whether
    each action happens depends on bits in X and Y which thus depends on the
    prior operations. However, we *can* determine, independently, whether each
    individual lane will swap X and Y and the XOR mask to apply to both before
    being passed on to the next lane pair. These values can be combined
    associatively and commutitavely, which makes it possible to determine the
    input transform for each lane pair in `log2(number_of_lane_pairs)` time.
    */
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
            let r_mask: i32x4 = r << (32 - $s) >> $s;
            d += ((3 * r.extract(0) as u64) ^ r.extract(1) as u64) << $s;

            let pos_inverted = n1 - pos;
            let newpos = bool32ix4::from_repr(r_mask.shuf(0, 0, 0, 0))
                .select(pos_inverted, pos).shuf(1, 0, 0, 0);
            pos = bool32ix4::from_repr(r_mask.shuf(1, 1, 1, 1))
                .select(pos, newpos);

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
}

pub fn hilbert_sort(array: &mut [CommonObject]) {
    array.sort_by_key(|o| xy_to_hilbert(o.p));
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
