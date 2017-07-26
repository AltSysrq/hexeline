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

/*! Facilities for doing 2D transforms and trigonometry.

The "fast sine" algorithm here is largely adapted from
https://web.archive.org/web/20080228213915/http://www.devmaster.net/forums/showthread.php?t=5784

Due to the uncertain availability of that page, here's a summary:

We want to approximate sin(θ) with a quadratic equation, i.e.,
```text
sin(θ) = Aθ²+Bθ+C
```

We know that sin(θ) = 0 and sin(π/2) = 1 and sin(π) = 0, which gives us the
three constraints needed to solve for the three variables:

```text
A0²+B0+C=0              C = 0
Aπ²/4+Bπ/2+C=1          Aπ²/4 + Bπ = 1
Aπ²+Bπ+C=0              Aπ²+Bπ = 0

Aπ² = -Bπ
Aπ = -B
A = -B/π

(-B/π)π²/4 + Bπ = 1
-Bπ/4 + Bπ/2 = 1
π/4 B = 1
B = 4/π

A = -4/π²
```

Thus, sin(θ) = -4/π² θ² + 4/π θ.

This works reasonably well for [0,π], but extends down to -∞ on [-π,0]. We
essentially want to mirror the parabola for negative θ. The second term gets
negated by nature of `θ` being negative, but we need to negate the first one
manually, which gives us a second equation, sin(θ) = 4/π² θ² + 4/π θ.

Unfortunately, this would mean having a branch in the code. However, we
actually get this for free if we evaluate θ² as θ|θ|, which equals θ² when θ>=0
but -θ² when θ<0. So our final quadratic formula is

```text
  sin(θ) = -4/π² θ|θ| + 4/π θ
```

This approximation is OK, but always overestimates sin(θ) except at the 3 fixed
points we originally set. We can get a better approximation by quadratically
adjusting the first approximation:

```text
  let y = first approximation
  sin(θ) = Qy + Py²
  Q + P = 1
```

The values of (Q,P) which produce the minimum absolute error are (0.775,0.225).
We also need to do the same absolute-value trick for this to work on both sides
of the domain. The formula can also be simplified:

```text
  sin(θ) = Py|y| + Qy
         = Py|y| + (1-P)y
         = Py|y| + y - Py
         = y(P|y| + 1 - P)
         = y(P|y| + Q)
```

(End of summary)

A lot of things are more convenient for us here since we use integers and
fixed-point. We don't need to worry about wrapping θ manually, because it
happens implicitly due to it being stored in a 16-bit integer. Since π = 32768,
the first approximation is mostly bit-shifts, with the only multiply being
θ|θ|. For the second approximation, we sacrifice a little bit of accuracy to
set P=1/4, which replaces more multiplies by bit-shifts.

This approach to cos/sin really is fast, too. On my main test system, computing
the full rotation matrix takes 3ns, whereas `f32::cos()` takes 7ns.

Besides basic trigonometry, another thing we need to do here is mimic cartesian
transformation in hexagonal space; i.e., compute a 2x2 matrix to produce the
same transform on (A,B) that as the original would in the (X,Y) coordinate
would have. We can find this by substituting a translation into cartesian
coordinates, rotation, and translation back. Note that the `k` term in the
transform (see `coords.rs`) always cancels out, so it is not notated here.

```text
  ; Translate new cartesian coords back to hexagonal
  a' = 2/3 * x'
  b' = -1/3 * x' + 1/sqrt(3) * y'

  | a' |   | 2/3        0          | | x' |
  | b' | = | -1/3       1/sqrt(3)  | | y' |

  ; Transform cartesian coords
  | x' |   | m0         m1         | | x |
  | y' | = | m2         m3         | | y |

  ; Transform from original hexagonal coords
  x = a - 1/2 * (b + c)
  y = sqrt(3)/2 * (b - c)
  c = -a - b
  x = a + 1/2 * a
  x = 3/2 * a
  y = sqrt(3)/2 * (2*b + a)
  y = sqrt(3)/2 * a + sqrt(3) * b

  | x |   | 3/2         0       | | a |
  | y | = | sqrt(3)/2   sqrt(3) | | b |

  Combine everything:

  | a' |   | 2/3   0         | | m0  m1 | | 3/2        0       | | a |
  | b' | = | -1/3  1/sqrt(3) | | m2  m3 | | sqrt(3)/2  sqrt(3) | | b |

  Premultiply the first two matrices:

  | a' |   | 2/3*m0                  2/3*m1                   | | 3/2        0       | | a |
  | b' | = | -1/3*m0 + 1/sqrt(3)*m2  -1/3*m1 + 1/sqrt(3) * m3 | | sqrt(3)/2  sqrt(3) | | b |

  Premultiply again:

  | a' |   | m0 + 1/sqrt(3)*m1                                 2/sqrt(3)*m1       | | a |
  | b' | = | -1/2*m0 + sqrt(3)/2*m2 - 1/sqrt(3)/2*m1 + 1/2*m3  -1/sqrt(3)*m1 + m3 | | b |
```

*/

use std::marker::PhantomData;
use std::ops;

use simd::*;
use simdext::*;

use physics::{Angle, DEG_90_CW, DEG_180};
use physics::coords::*;

/// A 2D affine transform.
///
/// This is a 2x2 matrix stored in row-major order. The fixed-point shift is
/// given by `AFFINE_POINT`.
#[derive(Debug, Clone, Copy)]
pub struct Affine2d<S : Space>(i32x4, PhantomData<S>);
pub type Affine2dO = Affine2d<Orthogonal>;
pub type Affine2dH = Affine2d<Hexagonal>;

/// Position of the fixed point in an `Affine2d`.
///
/// A value of 10 leaves 5 bits for the integer part.
pub const AFFINE_POINT: u32 = 10;

impl<S : Space> Default for Affine2d<S> {
    fn default() -> Affine2d<S> {
        Affine2d::identity()
    }
}

impl<S : Space> Affine2d<S> {
    #[inline]
    pub fn identity() -> Self {
        Affine2d(i32x4::new(1, 0, 0, 1), PhantomData)
    }

    #[inline(always)]
    pub fn from_repr(repr: i32x4) -> Self {
        Affine2d(repr, PhantomData)
    }

    #[inline(always)]
    pub fn repr(self) -> i32x4 {
        self.0
    }

    #[inline(always)]
    pub fn scale(x: i32, y: i32) -> Self {
        Affine2d(i32x4::new(x, 0, 0, y), PhantomData)
    }
}

impl Affine2dO {
    /// Computes the transform to rotate clockwise around the origin.
    pub fn rotate(sin_angle: Angle) -> Affine2dO {
        let cos_angle = sin_angle + DEG_90_CW;
        let neg_sin_angle = sin_angle + DEG_180;
        let theta = i32x4::new(cos_angle.0 as i32, neg_sin_angle.0 as i32,
                               sin_angle.0 as i32, cos_angle.0 as i32);
        // y = -4/π² θ|θ| + 4/π θ
        // π = 2**15
        // Y = -4*(2**POINT)/(2**15)² θ|θ| + 4*(2**POINT)/2**15 θ
        //   = -2**(POINT+2)/2**30 θ|θ| + 2**(POINT-13)θ
        //   = -2**(POINT-28) θ|θ| + 2**(POINT-13)θ
        //
        // In all calculations here, our fixed-point values are in the range
        // [-32768,+32768], so a 32-bit multiply won't overflow, and so we can
        // delay bit-shifting the point back to the correct place until after
        // the multiply.
        let first_term = theta * theta.abs() >> (28 - AFFINE_POINT);
        let second_term = theta >> (13 - AFFINE_POINT);
        let y = second_term - first_term;
        // Use the Py|y| + y - Py form so we don't need to load any other simd
        // constants. P=1/4, so Py = y >> 2
        let py = y >> 2;
        Affine2d(((py * y.abs()) >> AFFINE_POINT) + y - py, PhantomData)
    }

    /// Converts this transform so that it performs the same transform in
    /// hexagonal space.
    pub fn to_hexagonal(self) -> Affine2dH {
        /*
        | m0 + 1/sqrt(3)*m1                                 2/sqrt(3)*m1       |
        | -1/2*m0 + sqrt(3)/2*m2 - 1/sqrt(3)/2*m1 + 1/2*m3  -1/sqrt(3)*m1 + m3 |

        = m0            2/sqrt(3)*m1    sqrt(3)/2*m2    m3
        + 1/sqrt(3)*m1  0               -1/2*m0         -1/sqrt(3)*m1
        + 0             0               -1/sqrt(3)/2*m1 0
        + 0             0               1/2*m3          0

        m0 distinct factors: 1, -1/2
        m1 distinct factors: 1/sqrt(3), 2/sqrt(3), -1/sqrt(3)/2, -1/sqrt(3)
        m2 distinct factors: sqrt(3)/2
        m3 distinct factors: 1, 1/2

        Reorganise:

        = m0            0               1/2*m3          m3
        + 1/sqrt(3)*m1  2/sqrt(3)*m1    -1/sqrt(3)/2*m1 -1/sqrt(3)*m1
        + 0             0               sqrt(3)/2*m2   0
        - 0             0               1/2*m0          0

        Eliminate the -1/sqrt(3)/2*m1 term by moving it to a separate
        subtraction step:

        = m0            0               1/2*m3          m3
        + 1/sqrt(3)*m1  2/sqrt(3)*m1    sqrt(3)/2*m2   -1/sqrt(3)*m1
        - 0             0               1/sqrt(3)/2*m1  0
        - 0             0               1/2*m0          0

        Remove 0 from first row so `blend` can be used:

        = m0            1/2*m1          1/2*m3          m3
        + 1/sqrt(3)*m1  2/sqrt(3)*m1    sqrt(3)/2*m2   -1/sqrt(3)*m1
        - 0             0               1/sqrt(3)/2*m1  0
        - 0             1/2*m1          1/2*m0          0
        */
        let unit = self.0;
        let half: i32x4 = unit >> 1;
        let mulled: i32x4 = unit.shuf(1, 1, 2, 1) * i32x4::new(
            18919, 37837, 28378, -18919) >> 15;
        let half_mulled: i32x4 = mulled >> 1;

        Affine2d(i32x4::new(unit.extract(0), half.extract(1),
                            half.extract(3), unit.extract(3)) +
                 mulled -
                 i32x4::new(0, 0, half_mulled.extract(0), 0) -
                 i32x4::new(0, half.extract(1), half.extract(0), 0),
                 PhantomData)
    }
}

impl Affine2dH {
    /// Approximately equivalent to `Affine2dO::rotate(angle).to_hexagonal()`,
    /// but more efficient.
    pub fn rotate_hex(angle: Angle) -> Affine2dH {
        /*
        From to_hexagonal:

        = m0            0               1/2*m3          m3
        + 1/sqrt(3)*m1  2/sqrt(3)*m1    sqrt(3)/2*m2   -1/sqrt(3)*m1
        - 0             0               1/sqrt(3)/2*m1  0
        - 0             0               1/2*m0          0

        m0 = cos, m1 = -sin, m2 = sin, m3 = cos

        = cos           0               cos/2           cos
        + -sin/sqrt(3)  -2*sin/sqrt(3)  sin*sqrt(3)/2   sin/sqrt(3)
        + 0             0               sin/sqrt(3)/2   0
        + 0             0               -cos/2          0

        sin*sqrt(3)/2 + sin/sqrt(3)/2 =
        sin*3/sqrt(3)/2 + sin/sqrt(3)/2 =
        sin*4/sqrt(3)/2 =
        2*sin/sqrt(3)

        = cos           -sin/sqrt(3)    sin/sqrt(3)     cos
        + -sin/sqrt(3)  -sin/sqrt(3)    sin/sqrt(3)     sin/sqrt(3)

        Distinct factors:
        cos: 1
        sin: 1/sqrt(3)
        -sin: 1/sqrt(3)
        */
        let unit = Affine2dO::rotate(angle).repr();
        let mult: i32x4 = unit * i32x4::new(
            32768, 18919, 18919, 32768) >> 15;

        Affine2dH::from_repr(mult + mult.shuf(1, 1, 2, 2))
    }

    /// Assuming this transform was produced by `rotate_hex(theta)`, compute
    /// `rotate_hex(-theta)`.
    #[inline(always)]
    pub fn inv_rotate_hex(self) -> Self {
        /*
        cos(-θ) = cos(θ), sin(-θ) = -sin(θ)

        We're given

        [ cos-sin/sqrt(3)       -2sin/sqrt(3)
          2sin/sqrt(3)          cos+sin/sqrt(3) ]

        We want to get to

        [ cos+sin/sqrt(3)       2sin/sqrt(3)
          -2sin/sqrt(3)         cos-sin/sqrt(3) ]

        So we can simply permute the matrix.
        */
        Affine2dH::from_repr(self.repr().shuf(3, 2, 1, 0))
    }
}

impl<S : Space> ops::Mul<DualVector<S>> for Affine2d<S> {
    type Output = DualVector<S>;
    #[inline(always)]
    fn mul(self, v: DualVector<S>) -> DualVector<S> {
        let prod = v.repr().mulfp(self.0, AFFINE_POINT);
        DualVector::from_repr(prod.shuf(0, 2, 1, 3) +
                              prod.shuf(1, 3, 0, 2))
    }
}

impl<S : Space> ops::Mul<LineVector<S>> for Affine2d<S> {
    type Output = LineVector<S>;
    #[inline(always)]
    fn mul(self, v: LineVector<S>) -> LineVector<S> {
        let prod0202 = v.repr().shuf(0, 0, 2, 2).mulfp(
            self.0.shuf(0, 2, 0, 2), AFFINE_POINT);
        let prod1313 = v.repr().shuf(1, 1, 3, 3).mulfp(
            self.0.shuf(1, 3, 1, 3), AFFINE_POINT);
        LineVector::from_repr(prod0202 + prod1313)

    }
}

impl<S : Space> ops::Mul<Affine2d<S>> for Affine2d<S> {
    type Output = Affine2d<S>;

    #[inline(always)]
    fn mul(self, rhs: Affine2d<S>) -> Affine2d<S> {
        // | 00+12 01+13 |
        // | 20+32 21+33 |
        //
        // < 00 01 20 21 > +
        // < 12 13 32 33 >
        Affine2d(
            (self.0.shuf(0, 0, 2, 2) * rhs.0.shuf(0, 1, 0, 1) +
             self.0.shuf(1, 1, 3, 3) * rhs.0.shuf(2, 3, 2, 3)) >> AFFINE_POINT,
            PhantomData)
    }
}

#[cfg(test)]
mod test {
    use std::num::Wrapping;
    use test::{Bencher, black_box};

    use simdext::*;

    use physics::{DEG_90_CW, DEG_90_CCW, DEG_180};
    use super::*;

    #[test]
    fn simple_rotations() {
        let orig = Vod(1024, 256);
        let rot0 = Affine2d::rotate(Wrapping(0)) * orig;
        assert!(orig.repr().dist_2L1(rot0.repr()) < 8,
                "rotate(0) * {:?} => {:?}", orig, rot0);

        let rotcw = Affine2d::rotate(DEG_90_CW) * orig;
        assert!(Vod(-256, 1024).repr().dist_2L1(rotcw.repr()) < 8,
                "rotate(90) * {:?} => {:?}", orig, rotcw);

        let rotccw = Affine2d::rotate(DEG_90_CCW) * orig;
        assert!(Vod(256, -1024).repr().dist_2L1(rotccw.repr()) < 8,
                "rotate(-90) * {:?} => {:?}", orig, rotccw);

        let flip = Affine2d::rotate(DEG_180) * orig;
        assert!(Vod(-1024, -256).repr().dist_2L1(flip.repr()) < 8,
                "rotate(180) * {:?} => {:?}", orig, flip);
    }

    #[test]
    fn fine_rotation() {
        let orig = Vod(1024, 256);
        let rotated = Affine2d::rotate(Wrapping(5461 /* 30 deg */)) * orig;
        assert!(Vod(759, 734).repr().dist_2L1(rotated.repr()) < 16,
                "rotate(30) * {:?} => {:?}", orig, rotated);
    }

    #[test]
    fn matrix_mult() {
        let orig = Vod(1024, 256);
        let result = Affine2d::scale(2 << AFFINE_POINT, 1 << AFFINE_POINT) *
            Affine2d::rotate(DEG_90_CW) * orig;

        assert!(Vod(-512, 1024).repr().dist_2L1(result.repr()) < 8,
                "scale(2,1) * rotate(90) * {:?} => {:?}",
                orig, result);
    }

    #[test]
    fn line_mult() {
        let orig = Vol(1024, 2048, 4096, 8192);
        let xform = Affine2d::rotate(DEG_90_CW);

        let result = xform * orig;
        let expected_fst = xform * orig.fst_dual();
        let expected_snd = xform * orig.snd_dual();

        assert_eq!(result.x1(), expected_fst.x());
        assert_eq!(result.y1(), expected_fst.y());
        assert_eq!(result.x2(), expected_snd.x());
        assert_eq!(result.y2(), expected_snd.y());
    }

    #[test]
    fn affine_to_hexagonal() {
        let orig = Vod(102400, 25600);
        let rotation = Affine2d::rotate(Wrapping(5461 /* 30 deg */));
        let expected = rotation * orig;
        let actual = (rotation.to_hexagonal() * orig.single().to_vhr().dual())
            .to_vod();

        assert!(expected.repr().dist_2L1(actual.repr()) < 256,
                "Expected {:?}, got {:?}", expected, actual);
    }

    #[test]
    fn rotate_hex() {
        let orig = Vod(102400, 0);
        for theta in -32768i32..32768i32 {
            let theta = Wrapping(theta as i16);
            let rot_ortho = (Affine2d::rotate(theta) * orig).single()
                .to_vhr().single();
            let rot_hex = (Affine2d::rotate_hex(theta) *
                           orig.single().to_vhr().dual()).single();

            assert!(rot_ortho.repr().dist_2L1(rot_hex.repr()) < 256,
                    "{} rot({}) to_hex = {}, {} to_hex rot_hex({}) = {}, \
                     dist = {}", orig, theta.0, rot_ortho,
                    orig, theta.0, rot_hex,
                    rot_ortho.repr().dist_2L1(rot_hex.repr()));
        }
    }

    #[test]
    fn inv_rotate_hex() {
        for theta in -32768i32..32768i32 {
            let theta = Wrapping(theta as i16);
            let expected = Affine2dH::rotate_hex(-theta).repr();
            let actual = Affine2dH::rotate_hex(theta).inv_rotate_hex().repr();

            let diff = (expected - actual).abs();
            const SLOP: i32 = 4;
            assert!(diff.extract(0) <= SLOP &&
                    diff.extract(1) <= SLOP &&
                    diff.extract(2) <= SLOP &&
                    diff.extract(3) <= SLOP,
                    "For theta = {}, expected {:?}, got {:?}",
                    theta, expected, actual);
        }
    }

    #[bench]
    fn bench_rotate(b: &mut Bencher) {
        b.iter(|| Affine2d::rotate(black_box(Wrapping(1024))));
    }

    #[bench]
    fn bench_rotate_hex(b: &mut Bencher) {
        b.iter(|| Affine2d::rotate_hex(black_box(Wrapping(1024))));
    }

    #[bench]
    fn bench_matrix_vector_mult(b: &mut Bencher) {
        b.iter(|| black_box(Affine2dO::identity()) * black_box(Vod(1, 1)));
    }

    #[bench]
    fn bench_matrix_line_mult(b: &mut Bencher) {
        b.iter(|| black_box(Affine2dO::identity()) *
               black_box(Vol(1, 1, 2, 2)));
    }

    #[bench]
    fn bench_matrix_matrix_mult(b: &mut Bencher) {
        b.iter(|| black_box(Affine2dO::identity()) *
               black_box(Affine2dO::identity()));
    }

    // For comparison to our combined cos/sin
    #[bench]
    fn bench_libc_cos(b: &mut Bencher) {
        b.iter(|| black_box(1.2f32).cos())
    }
}
