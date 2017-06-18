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

/*! Utilities for working with coordinates.

## Coordinate Systems and Spaces

Internally, Hexeline uses no fewer than eight coordinate systems representing
two spaces.

The first space is _orthogonal space_, addressed with familiar cartesian
coordinates. X is to the right, Y is down. Z, for what it matters, is out of
the screen, so orthogonal space is left-handed. The Z coordinate is always 0.
Orthogonal space is generally only used for rendering or things outside
performance-sensitive areas where it is more convenient.

Everything performance-critical instead uses _hexagonal space_. When viewed as
a two-dimensional space, it is an _oblique_ space. As a three-dimensional
space, it is the plane X+Y+Z=0. For clarity, however, we refer to the
coordinates in hexagonal space as A, B, and C. A full description of hexagonal
space and why it is used is in a later section.

Within each space, there are four types of coordinates:

- Single coordinates ("S"). The first two dimensions are in SIMD lanes 0 and 1;
  lanes 2 and 3 have no values in particular.

- Redundant coordinates ("R"). All three dimensions are in the first three SIMD
  lanes, and the third SIMD lane has no particular value. This is "redundant"
  since the third coordinate can always be computed from the first two.

- Dual coordinates ("D"). The first coordinate is in SIMD lanes 0 *and* 2; the
  second is in lanes 1 and 3. This is the output of matrix-vector
  multiplication, and is more efficient to matrix-multiply.

- Line coordinates ("L"). The "coordinate" represents two points. The first two
  dimensions of the first point are in lanes 0 and 1; the first two dimensions
  of the second point are in lanes 2 and 3. Some operations can be applied to
  both points of L coordinates simultaneously.

The type for a vector in each space/coordinate pair is named with a `V`
followed by the space type (`o` or `h`) and then the coordinate type. For
example, `Vos` is orthogonal single coordinates; `Vhd` is hexagonal dual
coordinates.

## Hexagonal Space

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

use std::fmt;
use std::marker::PhantomData;
use std::mem;
use std::ops;

use simd::*;
use simdext::*;

/// The amount to left-shift hexagonal coordinates 1 by to get `CELL_RADIUS`.
pub const CELL_RADIUS_SHIFT: u8 = 8;
/// The  distance from the centre of a cell to any of its vertices.
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

pub trait Space {
    fn coord_prefix() -> char;
    fn coord_suffix() -> char;
    fn compute_z(x: i32, y: i32) -> i32;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Orthogonal;
impl Space for Orthogonal {
    fn coord_prefix() -> char { '[' }
    fn coord_suffix() -> char { ']' }
    #[inline(always)]
    fn compute_z(_x: i32, _y: i32) -> i32 { 0 }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Hexagonal;
impl Space for Hexagonal {
    fn coord_prefix() -> char { '<' }
    fn coord_suffix() -> char { '>' }
    #[inline(always)]
    fn compute_z(a: i32, b: i32) -> i32 { -a-b }
}


pub struct SingleVector<S : Space>(i32x4, PhantomData<S>);
pub struct RedundantVector<S : Space>(i32x4, PhantomData<S>);
pub struct DualVector<S : Space>(i32x4, PhantomData<S>);
pub struct LineVector<S : Space>(i32x4, PhantomData<S>);

impl<S : Space> SingleVector<S> {
    #[inline(always)]
    pub fn new(x: i32, y: i32) -> Self {
        SingleVector(i32x4::new(x, y, unsafe { mem::uninitialized() },
                                unsafe { mem::uninitialized() }),
                     PhantomData)
    }

    #[inline(always)]
    pub fn redundant(self) -> RedundantVector<S> {
        RedundantVector::new(
            self.0.extract(0), self.0.extract(1),
            S::compute_z(self.0.extract(0), self.0.extract(1)))
    }

    #[inline(always)]
    pub fn line(self, two: Self) -> LineVector<S> {
        LineVector::new(self.0.extract(0), self.0.extract(1),
                        two.0.extract(0), two.0.extract(1))
    }

    #[inline(always)]
    pub fn dual(self) -> DualVector<S> {
        DualVector::new(self.0.extract(0), self.0.extract(1))
    }
}

impl<S : Space> RedundantVector<S> {
    #[inline(always)]
    pub fn new(x: i32, y: i32, z: i32) -> Self {
        RedundantVector(i32x4::new(x, y, z, unsafe { mem::uninitialized() }),
                        PhantomData)
    }

    #[inline(always)]
    pub fn single(self) -> SingleVector<S> {
        SingleVector(self.0, PhantomData)
    }

    #[inline(always)]
    pub fn dual(self) -> DualVector<S> {
        self.single().dual()
    }
}

impl<S : Space> DualVector<S> {
    #[inline(always)]
    pub fn new(x: i32, y: i32) -> Self {
        DualVector(i32x4::new(x, y, x, y), PhantomData)
    }

    #[inline(always)]
    pub fn single(self) -> SingleVector<S> {
        SingleVector(self.0, PhantomData)
    }

    #[inline(always)]
    pub fn redundant(self) -> RedundantVector<S> {
        self.single().redundant()
    }
}

impl<S : Space> LineVector<S> {
    #[inline(always)]
    pub fn new(x1: i32, y1: i32, x2: i32, y2: i32) -> Self {
        LineVector(i32x4::new(x1, y1, x2, y2), PhantomData)
    }

    #[inline(always)]
    pub fn fst(self) -> SingleVector<S> {
        SingleVector(self.0, PhantomData)
    }

    #[inline(always)]
    pub fn fst_dual(self) -> DualVector<S> {
        DualVector(self.0.shuf(0, 1, 0, 1), PhantomData)
    }

    #[inline(always)]
    pub fn snd(self) -> SingleVector<S> {
        SingleVector(self.0.shuf(2, 3, 2, 3), PhantomData)
    }

    #[inline(always)]
    pub fn snd_dual(self) -> DualVector<S> {
        DualVector(self.0.shuf(2, 3, 2, 3), PhantomData)
    }
}

macro_rules! binop {
    ($S:ident, $lhs:ident, $rhs:ty, $rhsn:ident, $rhs_access:expr,
     $name:ident, $meth:ident, $op:tt) => {
        impl<$S : Space> ops::$name<$rhs> for $lhs<$S> {
            type Output = Self;
            #[inline(always)]
            fn $meth(self, $rhsn: $rhs) -> Self {
                $lhs(self.0 $op $rhs_access, PhantomData)
            }
        }
    }
}

macro_rules! vector_common {
    ($name:ident $(,$ocoord:ident / $hcoord:ident = $ix:tt)+) => {
        impl<S : Space> Clone for $name<S> {
            #[inline(always)]
            fn clone(&self) -> Self {
                $name(self.0, PhantomData)
            }
        }

        impl<S : Space> Copy for $name<S> { }

        impl<S : Space> $name<S> {
            #[inline(always)]
            pub fn repr(self) -> i32x4 {
                self.0
            }

            #[inline(always)]
            pub fn from_repr(repr: i32x4) -> Self {
                $name(repr, PhantomData)
            }

            #[inline(always)]
            pub fn multp(self, fp: i32x4, point: u32) -> Self {
                $name(self.0.mulfp(fp, point), PhantomData)
            }
        }

        impl $name<Orthogonal> {
            $(#[inline(always)]
              pub fn $ocoord(self) -> i32 {
                  self.0.extract($ix)
              })*
        }
        impl $name<Hexagonal> {
            $(#[inline(always)]
              pub fn $hcoord(self) -> i32 {
                  self.0.extract($ix)
              }
            )*
        }

        impl<S : Space> fmt::Debug for $name<S> {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                let _ch = S::coord_prefix();
                $(write!(f, "{}{}", _ch, self.0.extract($ix))?; let _ch = ',';)*
                write!(f, "{}", S::coord_suffix())
            }
        }

        impl<S : Space> fmt::Display for $name<S> {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                fmt::Debug::fmt(self, f)
            }
        }

        binop!(S, $name, $name<S>, rhs, rhs.0, Add, add, +);
        binop!(S, $name, $name<S>, rhs, rhs.0, BitAnd, bitand, &);
        binop!(S, $name, $name<S>, rhs, rhs.0, BitOr, bitor, |);
        binop!(S, $name, $name<S>, rhs, rhs.0, BitXor, bitxor, ^);
        binop!(S, $name, $name<S>, rhs, rhs.0, Mul, mul, *);
        binop!(S, $name, u8, rhs, rhs, Shl, shl, <<);
        binop!(S, $name, u32, rhs, rhs, Shl, shl, <<);
        binop!(S, $name, u8, rhs, rhs, Shr, shr, >>);
        binop!(S, $name, u32, rhs, rhs, Shr, shr, >>);
        binop!(S, $name, $name<S>, rhs, rhs.0, Sub, sub, -);

        impl<S : Space> ops::Neg for $name<S> {
            type Output = Self;
            #[inline(always)]
            fn neg(self) -> Self {
                $name(i32x4::splat(0) - self.0, PhantomData)
            }
        }

        impl<S : Space> Default for $name<S> {
            #[inline(always)]
            fn default() -> Self {
                $name(i32x4::splat(0), PhantomData)
            }
        }
    }
}

vector_common!(SingleVector, x/a = 0, y/b = 1);
vector_common!(RedundantVector, x/a = 0, y/b = 1, z/c = 2);
vector_common!(DualVector, x/a = 0, y/b = 1);
vector_common!(LineVector, x1/a1 = 0, y1/b1 = 1, x2/a2 = 2, y2/b2 = 3);

macro_rules! vtype {
    ($(pub type $name:ident = $class:ident<$space:ident>(
        $($coord:ident),*);)*) =>
    {$(
        pub type $name = $class<$space>;
        #[inline(always)]
        #[allow(non_snake_case)]
        pub fn $name($($coord: i32),*) -> $name {
            $class::new($($coord),*)
        }
    )*}
}
vtype! {
    pub type Vos = SingleVector<Orthogonal>(x, y);
    pub type Vor = RedundantVector<Orthogonal>(x, y, z);
    pub type Vod = DualVector<Orthogonal>(x, y);
    pub type Vol = LineVector<Orthogonal>(x1, y1, x2, y2);
    pub type Vhs = SingleVector<Hexagonal>(a, b);
    pub type Vhr = RedundantVector<Hexagonal>(a, b, c);
    pub type Vhd = DualVector<Hexagonal>(a, b);
    pub type Vhl = LineVector<Hexagonal>(a1, b1, a2, b2);
}

impl Vos {
    pub fn to_vhr(self) -> Vhr {
        let x = i32x4::splat(self.x());
        let y = i32x4::splat(self.y());
        let col1 = x.mulfp(i32x4::new(26755, -13377, -13377, 0), 15);
        let col2 = y.mulfp(i32x4::new(0, 23170, -23170, 0), 15);
        Vhr::from_repr(col1 + col2)
    }
}

impl Vhs {
    /// Return the (A,B) index of the hexagonal cell containing this
    /// coordinate.
    #[inline(never)]
    pub fn to_index(self) -> (i32, i32) {
        // The core idea here is that we'll be rounding A and B up or down to the
        // nearest integer. They don't necessarily round the same way, so there are
        // 4 options total. C gets rounded similarly, though we consider it in
        // terms of A and B instead. Each rounding gives us a signed "residue";
        // correct solutions are those where the maximum residue among the three
        // axes minus the minimum residue of same is less than or equal to the cell
        // diameter.

        // Put A and B in their own registers
        let base_a = i32x4::splat(self.a());
        let base_b = i32x4::splat(self.b());
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

        let max_residue = a_residue.nsw_max(b_residue).nsw_max(c_residue);
        let min_residue = a_residue.nsw_min(b_residue).nsw_min(c_residue);
        // We can't use .le() because it gets turned into a function call even with
        // SSE4.1. We can accomplish (a <= b) with ((a - b - 1) >> 31) though
        // (which also gets us a convenient bitmask). Though here we only care
        // about the sign bit, so the >> 31 is elided.
        let valid_mask: i32x4 =
            (max_residue - min_residue) - i32x4::splat(CELL_COORD_MASK + 2);
        // Efficiently get a 4-bit value indicating which lanes of valid_mask have
        // their sign bit set (i.e., are valid solutions).
        let valid_bits = valid_mask.movemask();
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

        let coord = self.repr() >> CONTINUOUS_TO_CELL_SHIFT;
        (coord.extract(0) + a_off, coord.extract(1) + b_off)
    }
}

impl Vhd {
    pub fn to_vod(self) -> Vod {
        /*
        x = k*a - k/2*(c + b)
        y = k*sqrt(3)/2 * (b - c)

        c = -a-b
        x = k*a - k/2*(-a)
        x = 3/2*k*a
        y = k*sqrt(3)/2 * (a + 2*b)
        y = k*sqrt(3)/2 * a + k*sqrt(3) * b

        k = sqrt(2)/sqrt(3)
        x = sqrt(3)/sqrt(2) * a
        y = 1/sqrt(2) * a + sqrt(2) * b

        | x |   | 3/2        0       |
        | y | = | 1/sqrt(2)  sqrt(2) |
        */
        let prod = self.repr().mulfp(
            i32x4::new(20066, 0, 11585, 23170), 14);
        Vod::from_repr(prod.shuf(0, 2, 0, 2) + prod.shuf(1, 3, 1, 3))
    }
}

#[cfg(test)]
mod test {
    use std::hash::Hasher;
    use test::{Bencher, black_box};

    use fnv;
    use simdext::*;

    use super::*;

    macro_rules! assert_approx {
        ($fuzz:expr, $expected:expr, $actual:expr) => {{
            let expected = $expected;
            let actual = $actual;
            assert!(actual >= expected - $fuzz && actual <= expected + $fuzz,
                    "Expected {:?} +/- {:?}, got {:?}",
                    expected, $fuzz, actual);
        }}
    }

    #[bench]
    fn bench_vos_to_vhr(b: &mut Bencher) {
        b.iter(|| black_box(Vos(65536, 65536)).to_vhr())
    }

    #[bench]
    fn bench_vhs_to_index(b: &mut Bencher) {
        b.iter(|| black_box(Vhs(65536, 65536)).to_index())
    }

    #[test]
    fn vos_to_vhr_smoke_test() {
        // Check a few basic properties. This doesn't exhaustively check that
        // the code is correct. The implementation here was verified by
        // graphing the results and manually verifying them.
        let zero = Vos(0, 0).to_vhr();
        assert_eq!(0, zero.a());
        assert_eq!(0, zero.b());
        assert_eq!(0, zero.c());

        let pos_x = Vos(256, 0).to_vhr();
        assert!(pos_x.a() > 128, "Expected {} > 128", pos_x.a());
        assert!(pos_x.b() < 0);
        assert!(pos_x.c() < 0);
        assert_approx!(256, 65536, pos_x.repr().dist_3L2_squared(zero.repr()));
        assert_approx!(8, 0, pos_x.repr().hsum_3());

        let pos_y = Vos(0, 256).to_vhr();
        assert!(pos_y.b() > 0);
        // Expect 5 bits of precision to be lost
        assert_approx!(512, 65536, pos_y.repr().dist_3L2_squared(zero.repr()));
        assert_approx!(8, 0, pos_y.repr().hsum_3());
    }

    #[test]
    fn vos_to_vhr_reproducibility() {
        // Ensure that all implementations (i.e., compilations with different
        // CPU features) produce the exact same results.
        let mut hasher = fnv::FnvHasher::default();
        for y in -1000..1000 {
            for x in -1000..1000 {
                let hexa = Vos(x, y).to_vhr();
                hasher.write_i32(hexa.a());
                hasher.write_i32(hexa.b());
            }
        }

        assert_eq!(8664580250012815522, hasher.finish());
    }

    #[test]
    fn vhs_to_index_smoke_test() {
        let zero = Vhs(0, 0).to_index();
        assert_eq!(0, zero.0);
        assert_eq!(0, zero.1);

        let pos_x = Vos(16384, 0).to_vhr().single().to_index();
        assert!(pos_x.0 > 16, "Expected A = {} to be > 16", pos_x.0);
        assert!(pos_x.1 < -8, "Expected B = {} to be < -8", pos_x.1);
    }

    #[test]
    fn hexagonal_to_index_reproducibility() {
        let mut hasher = fnv::FnvHasher::default();
        for y in -1000..1000 {
            for x in -1000..1000 {
                let index = Vhs(x, y).to_index();
                hasher.write_i32(index.0);
                hasher.write_i32(index.1);
            }
        }

        assert_eq!(5086880218595302901, hasher.finish());
    }

    #[test]
    fn vhd_to_vod_inverse() {
        // Test hexagonal_to_cartesian in terms of inverting
        // cartesian_to_hexagonal and also ensure reproduciblity.
        let mut hasher = fnv::FnvHasher::default();

        for y in -1000..1000 {
            for x in -1000..1000 {
                let cart = Vos(x, y);
                let hexa = cart.to_vhr();
                let cart2 = hexa.dual().to_vod();
                assert!(cart.repr().dist_2L1(cart2.repr()) <= 32,
                        "{} => {} => {} (dist {})",
                        cart, hexa, cart2,
                        cart.repr().dist_2L1(cart2.repr()));

                hasher.write_i32(cart2.x());
                hasher.write_i32(cart2.y());
            }
        }

        assert_eq!(12864379608347279471, hasher.finish());
    }
}
