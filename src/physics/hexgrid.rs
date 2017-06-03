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

/*! Utilities for working in hexagonal space.

Inspiration for the core idea here comes from
http://keekerdc.com/2011/03/hexagon-grids-coordinate-systems-and-distance-calculations/

The root of the problem we are trying to solve here is to be able to address
positions in composites (i.e., hexagonal grids) via continuous external 2D
cartesian coordinates.

The core idea is to translate (X,Y) coordinates into 3-space. For the latter,
we use the coordinate names (A,B,C) for clarity. We constrain 3-space to the
plane `A + B + C = 0` and project it into 2-space with an orthographic position
such that each 3-space axis is at a 120° angle to the other two. The key here
is that this allows to align the axes with the edges of the hexagons. We thus
get something like this:

```text
*-----> +X
|   ^ +C
|    \
|     \
|      \
|       \
v        *------------> +A
+Y      / ___________      \
       / /     |     \     |
      / /      |      \    | CELL_H
     / /       |       \   |
    v /_________________\  /
   +B \        |        /
       \       |       /
        \      |      /
      |  \_____|_____/
      |
      |        \____/ CELL_CORNER
      \________/

      CELL_RADIUS
```

Every cell is indexed by an (A,B,C) coordinate where all axes are integer
multiples of `2*CELL_RADIUS`. This coordinate is the centre of the hexagon.
Since C can be trivially derived from A and B, we can store hexagons in a
two-dimensional array addressed directly by the (A,B) pair.

Transforming a 2-space point to hex-space is a simple affine transform. We
derive it here long-hand. `k` is a scaling factor we'll address later.

  Given:
  x = k * a - k/2 * (c + b)
  y = k*sqrt(3)/2 * (b - c)
  a + b + c = 0

  Work:
  b = -a - c
  x = k * a - k/2 * (c - a - c) = k * a + k/2 * a = 3*k/2 * a
  a = 2/3/k * x

  b = -2/3/k * x - c
  y = k*sqrt(3)/2 * (-2/3/k * x - 2*c)
  -2*y/k/sqrt(3) = 2/3/k * x + 2*c
  2*c = -2*y/k/sqrt(3) - 2/3/k * x
  c = -y/k/sqrt(3) - 1/3/k * x

  b = -2/3/k * x + y/k/sqrt(3) + 1/3/k * x
  b = -1/3/k * x + y/k/sqrt(3)

  Final transformation:

    | a |   | 2/3/k     0               c1 | | x |
    | b | = | -1/3/k    1/k/sqrt(3)     c2 | | y |
    | c |   | -1/3/k    -1/k/sqrt(3)    c3 | | 0 |

Three of the values in the matrix can have any arbitrary value since the Z
coordinate of the input is always 0.

We want the determinate of the transformation to be 1 so that the L2 metric in
hexagonal space equals the L2 metric of the equivalent coordinates in cartesian
space. However, we the determinate depends on c3, which is also unconstrained,
so we can't solve it this way.

Instead, go directly from the distance equations:

  Goal: sqrt(a**2 + b**2 + c**2) = sqrt(x**2 + y**2)
        a**2 + b**2 + c**2 = x**2 + y**2

  x**2 + y**2
  (k * a - k/2 * (c + b)) ** 2 + (k*sqrt(3)/2 * (b - c)) ** 2

  k**2 * ((a - 1/2*(c+b))**2 + 3/4 * (b-c)**2)
  k**2 * ((a**2 - a*(c+b) + 1/4*(c+b)**2) +
          (3/4 * b**2 - 3/2 * b * c + 3/4 * c**2))

  k**2 * (a**2 - a*b - a*c + 1/4*c**2 + 1/2*b*c + 1/4*b**2 +
          3/4 * b**2 - 3/2*b*c + 3/4*c**2)

  k**2 * (a**2 + b**2 + c**2 - a*b - a*c - b*c)
  k**2 * (a**2 + b**2 + c**2 - (a*b + a*c + b*c))

  Looking at that last group only:
  a*b + a*c + b*c
  Since a = -(b+c), b=-(a+c), c=-(a+b):
  -b*(b+c) - a*(a+b) - c*(a+c)
  -b**2 - b*c - a**2 - a*b - c**2 - a*c

  Thus
  a*b + a*c + b*c = -a**2 - b**2 - c**2 - a*b - a*c - b*c
  2*(a*b + a*c + b*c) = -a**2 - b**2 - c**2
  a*b + a*c + b*c = -1/2*a**2 -1/2*b**2 - 1/2*c**2

  Substituting back in:
  k**2 * (a**2 + b**2 + c**2 - (-1/2*a**2 -1/2*b**2 - 1/2*c**2))
  k**2 * (3/2 * a**2 + 3/2 * b**2 + 3/2 * c**2)
  k**2 * 3/2 * (a**2 + b**2 + c**2)

  So
  k**2 * 3/2 * (a**2 + b**2 + c**2) = a**2 + b**2 + c**2
  k**2 * 3/2 = 1
  k**2 = 2/3
  k = sqrt(2)/sqrt(3)

Thus the substituted transformation matrix is:
    | a | = | 2*sqrt(3)/sqrt(2)/3  0          c1 | | x |
    | b | = | -sqrt(3)/sqrt(2)/3   1/sqrt(2)  c2 | | y |
    | c | = | -sqrt(3)/sqrt(2)/3  -1/sqrt(2)  c3 | | 0 |

It would also be nice for relative angles to be preserved. We get this for free
since this is a length-preserving transformation; in order to skew angles, it
would necessarily need to compress or expand space in some cases, but we know
from above that this never happens. I.e., if we have an arbitrary triangle and
subject it to the transform, the lengths of the edges of the transformed
triangle match the original triangle, thus the new angles of the triangle must
match as well.

A common thing to do is to translate from continuous hexagonal coordinates to
the containing hexagon. How to do this is not obvious.

The coordinate system can be counter-intuitive. Particularly, the (A,B)
coordinates can actually be _less_ than the hexagon's origin; for example, as
one moves straight along positive X, the B and C coordinates become _smaller_.
It helps to first visualise the local neighbourhood of a single cell:

```text
                   ____
                  /    \
             ____/ 0-11 \____
            /    \      /    \
           / -101 \____/ 1-10 \
           \      /    \      /
            \____/ 000  \____/
            /    \      /    \
           / -110 \____/ 10-1 \
           \      /    \      /
            \____/ 01-1 \____/
                 \      /
                  \____/

```

After a bunch of drawing hexagons on paper and staring at the numbers, you
eventually find that the boundary condition for each edge is

  max(dA,dB,dC) - min(dA,dB,dC) = 1

where `dA` is the A distance from the origin of the hexagon to the point, etc.

Substituting some stuff gives us

  max(dA,dB,-(dA+dB)) - min(dA,dB,-(dA+dB)) <= 1

as the interior of the hexagon.
*/

use std::mem;

use simd::*;

/// The amount to left-shift 1 by to get `CELL_RADIUS`.
pub const CELL_RADIUS_SHIFT: u8 = 8;
/// The distance from the centre of a cell to any of its vertices.
///
/// An outer radius of 512 means that the outer diameter is 1024, or 1/64th of
/// a screen width.
pub const CELL_RADIUS: i32 = 1 << CELL_RADIUS_SHIFT;
/// The distance from the centre of a cell to the centre of one of the edges.
///
pub const CELL_H: i32 = CELL_RADIUS * 866_026 / 1_000_000;
/// The distance from the centre of a cell edge to the vertex on either side.
pub const CELL_CORNER: i32 = CELL_RADIUS / 2;

/// The amount to right-shift hexagonal coordinates to go from continuous
/// coordinate to containing cell coordinate.
pub const CONTINUOUS_TO_CELL_SHIFT: u8 = CELL_RADIUS_SHIFT + 1;
/// Bitmask to position continuous hexagonal coordinates over the origin of the
/// containing cell.
pub const CELL_COORD_MASK: i32 = (1 << CONTINUOUS_TO_CELL_SHIFT) - 1;

/// Convert the given cartesian coordinates (X in component 0, Y in component
/// 1) to hexagonal coordinates (A, B, C in 0, 1, 2). Output component 3 is
/// always 0.
pub fn cartesian_to_hexagonal(cart: i32x4) -> i32x4 {
    let x = i32x4::splat(cart.extract(0));
    let y = i32x4::splat(cart.extract(1));
    let col1 = x * i32x4::new(6689, -3344, -3344, 0);
    let col2 = y * i32x4::new(0, 5792, -5792, 0);
    (col1 + col2) >> 13
}

#[inline(never)]
pub fn hexagonal_to_index(hexa: i32x4) -> (i32, i32) {
    // The core idea here is that we'll be rounding A and B up or down to the
    // nearest integer. They don't necessarily round the same way, so there are
    // 4 options total. C gets rounded similarly, though we consider it in
    // terms of A and B instead. Each rounding gives us a signed "residue";
    // correct solutions are those where the maximum residue among the three
    // axes minus the minimum residue of same is less than or equal to the cell
    // diameter.

    // Put A and B in their own registers
    let base_a = i32x4::splat(hexa.extract(0));
    let base_b = i32x4::splat(hexa.extract(1));
    // Calculate our rounding options. Note that for A we round up in the last
    // two lanes, where for B we round up in the odd numbered lanes.
    let rounded_a =
        (base_a + i32x4::new(0, 0, CELL_COORD_MASK+1, CELL_COORD_MASK+1)) &
        i32x4::splat(!CELL_COORD_MASK);
    let rounded_b =
        (base_b + i32x4::new(0, CELL_COORD_MASK+1, 0, CELL_COORD_MASK+1)) &
        i32x4::splat(!CELL_COORD_MASK);
    // Compute the residue for all for options, as well as for C
    let a_residue = rounded_a - base_a;
    let b_residue = rounded_b - base_b;
    // i32x4::neg expands to a function call, so subtract from 0 instead.
    let c_residue = i32x4::splat(0) - a_residue - b_residue;

    // TODO Support systems without SSE4.1
    use simd::x86::sse4_1::Sse41I32x4;
    let max_residue = a_residue.max(b_residue).max(c_residue);
    let min_residue = a_residue.min(b_residue).min(c_residue);
    // We can't use .le() because it gets turned into a function call even with
    // SSE4.1. We can accomplish (a <= b) with ((a - b - 1) >> 31) though
    // (which also gets us a convenient bitmask). Though here we only care
    // about the sign bit, so the >> 31 is elided.
    let valid_mask: i32x4 =
        (max_residue - min_residue) - i32x4::splat(CELL_COORD_MASK + 2);
    // Efficiently get a 4-bit value indicating which lanes of valid_mask have
    // their sign bit set (i.e., are valid solutions).
    // TODO Support non-SSE
    extern "platform-intrinsic"{
        fn x86_mm_movemask_ps(a: f32x4) -> i32;
    }
    let valid_bits = unsafe {
        x86_mm_movemask_ps(mem::transmute(valid_mask))
    } as u32;
    // Lookup-table-in-a-register
    // Map movemask values to A and B offsets.
    // We know there is always at least 1 solution, so if we bit-shift the
    // movemask output right one, we know that a value of 0 indicates that only
    // lane 0 was set, and so we can cram each of these into 8 bits.
    //                        7  6  5  4  3  2  1  0
    //                     3  1  1  1  1  0  0  0  ?
    //                     2  1  1  0  0  1  1  0  ?
    //                     1  1  0  1  0  1  0  1  ?
    //                     0  ?  ?  ?  ?  ?  ?  ?  1
    //                choose  3  3  3  3  2  2  1  0
    let valid_shuf_tab_a = 0b_1__1__1__1__1__1__0__0u16;
    let valid_shuf_tab_b = 0b_1__1__1__1__0__0__1__0u16;
    let shuf_ix = valid_bits >> 1;
    let a_off = (valid_shuf_tab_a >> shuf_ix) as i32 & 1;
    let b_off = (valid_shuf_tab_b >> shuf_ix) as i32 & 1;

    let coord = hexa >> CONTINUOUS_TO_CELL_SHIFT;
    (coord.extract(0) + a_off, coord.extract(1) + b_off)
}

/// Approximately inverts `cartesian_to_hexagonal`.
///
/// X and Y are stored in output components 0 and 1. Components 2 and 3 are zero.
pub fn hexagonal_to_cartesian(hexa: i32x4) -> i32x4 {
    /*
    x = k*a - k/2*(c + b)
    y = k*sqrt(3)/2 * (b - c)

    x = k*a - k/2*b - k/2*c
    y = k*sqrt(3)/2 * b - k*sqrt(3)/2 * c

    x = sqrt(2)/sqrt(3)*a - sqrt(2)/sqrt(3)*b - sqrt(2)/sqrt(3)*b
    y = sqrt(2)/2*b - sqrt(2)/2*c

    | x |   | sqrt(2)/sqrt(3)   -sqrt(2)/sqrt(3)        -sqrt(2)/sqrt(3) | | a |
    | y | = | 0                 sqrt(2)/2               -sqrt(2)/2       | | b |
    | _ |   | 0                 0                       0                | | c |
    */
    // Make a vector of (b,b,a,a) so we can do the first two columns at the
    // same time.
    let bbaa = i32x4::new(hexa.extract(1), hexa.extract(1),
                          hexa.extract(0), hexa.extract(0));
    // Compute the first two columns at the same time. The second row of the
    // first column (3 component) is simply set to the negative of the first
    // row, which we use to cancel it out later.
    let col21: i32x4 = (bbaa * i32x4::new(-6689, 5793, 6689, -6689)) >> 13;
    let c = i32x4::splat(hexa.extract(2));
    let col3 = (c * i32x4::new(-6689, 5793, 0, 0)) >> 13;

    col21 + col3 + i32x4::new(col21.extract(2), col21.extract(3),
                              // Cancel out remnants of the 2nd half of col21
                              col21.extract(3), col21.extract(2))
}

#[cfg(test)]
mod test {
    use super::*;

    use test::{Bencher, black_box};

    #[bench]
    fn bench_cartesian_to_hexagonal(b: &mut Bencher) {
        b.iter(|| cartesian_to_hexagonal(black_box(
            i32x4::splat(65536))))
    }

    #[bench]
    fn bench_hexagonal_to_index(b: &mut Bencher) {
        b.iter(|| hexagonal_to_index(black_box(
            i32x4::splat(65536))))
    }
}
