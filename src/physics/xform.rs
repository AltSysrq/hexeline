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

The values of (Q,P) which produce the minimum absolute error are (0.755,0.225).
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

Besides basic trigonometry, another thing we need to do here is mimic cartesian
transformation in hexagonal space; i.e., compute a 2x2 matrix to produce the
same transform on (A,B) that as the original would in the (X,Y) coordinate
would have. We can find this by substituting a translation into cartesian
coordinates, rotation, and translation back. Note that the `k` term in the
transform (see `hexgrid.rs`) always cancels out, so it is not notated here.

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

use std::ops;

use simd::*;
use simdext::*;

use physics::{Angle, DEG_90_CW, DEG_180};

/// A 2D affine transform.
///
/// This is a 2x2 matrix stored in row-major order. The fixed-point shift is
/// given by `AFFINE_POINT`.
#[derive(Debug, Clone, Copy)]
pub struct Affine2d(pub i32x4);

/// Position of the fixed point in an `Affine2d`.
///
/// A value of 10 leaves 5 bits for the integer part.
pub const AFFINE_POINT: u32 = 10;

impl Default for Affine2d {
    fn default() -> Affine2d {
        Affine2d::identity()
    }
}

impl Affine2d {
    pub fn identity() -> Affine2d {
        Affine2d(i32x4::new(1, 0, 0, 1))
    }

    pub fn rotate(sin_angle: Angle) -> Affine2d {
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
        Affine2d(((py * y.abs()) >> AFFINE_POINT) + y - py)
    }

    #[inline(always)]
    pub fn scale(x: i32, y: i32) -> Self {
        Affine2d(i32x4::new(x, 0, 0, y))
    }

    pub fn to_hexagonal(self) -> Self {
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
                 i32x4::new(0, half.extract(1), half.extract(0), 0))
    }
}

impl ops::Mul<i32x4> for Affine2d {
    type Output = i32x4;
    fn mul(self, v: i32x4) -> i32x4 {
        // TODO We'd need fewer shuffles here if things were laid out better
        let prod = v.shuf(0, 0, 1, 1).mulfp(
            self.0.shuf(0, 2, 1, 3), AFFINE_POINT);
        prod + prod.shuf(2, 3, 2, 3)
    }
}

impl ops::Mul<Affine2d> for Affine2d {
    type Output = Affine2d;

    fn mul(self, rhs: Affine2d) -> Affine2d {
        // | 00+12 01+13 |
        // | 20+32 21+33 |
        //
        // < 00 01 20 21 > +
        // < 12 13 32 33 >
        Affine2d(
            (self.0.shuf(0, 0, 2, 2) * rhs.0.shuf(0, 1, 0, 1) +
             self.0.shuf(1, 1, 3, 3) * rhs.0.shuf(2, 3, 2, 3)) >> AFFINE_POINT)
    }
}

#[cfg(test)]
mod test {
    use std::num::Wrapping;
    use test::{Bencher, black_box};

    use simd::*;
    use simdext::*;

    use physics::{DEG_90_CW, DEG_90_CCW, DEG_180};
    use physics::hexgrid;
    use super::*;

    #[test]
    fn simple_rotations() {
        let orig = i32x4::new(1024, 256, 0, 0);
        let rot0 = Affine2d::rotate(Wrapping(0)) * orig;
        assert!(orig.dist_2L1(rot0) < 8,
                "rotate(0) * {:?} => {:?}", orig, rot0);

        let rotcw = Affine2d::rotate(DEG_90_CW) * orig;
        assert!(i32x4::new(-256, 1024, 0, 0).dist_2L1(rotcw) < 8,
                "rotate(90) * {:?} => {:?}", orig, rotcw);

        let rotccw = Affine2d::rotate(DEG_90_CCW) * orig;
        assert!(i32x4::new(256, -1024, 0, 0).dist_2L1(rotccw) < 8,
                "rotate(-90) * {:?} => {:?}", orig, rotccw);

        let flip = Affine2d::rotate(DEG_180) * orig;
        assert!(i32x4::new(-1024, -256, 0, 0).dist_2L1(flip) < 8,
                "rotate(180) * {:?} => {:?}", orig, flip);
    }

    #[test]
    fn fine_rotation() {
        let orig = i32x4::new(1024, 256, 0, 0);
        let rotated = Affine2d::rotate(Wrapping(5461 /* 30 deg */)) * orig;
        assert!(i32x4::new(759, 734, 0, 0).dist_2L1(rotated) < 16,
                "rotate(30) * {:?} => {:?}", orig, rotated);
    }

    #[test]
    fn matrix_mult() {
        let orig = i32x4::new(1024, 256, 0, 0);
        let result = Affine2d::scale(2 << AFFINE_POINT, 1 << AFFINE_POINT) *
            Affine2d::rotate(DEG_90_CW) * orig;

        assert!(i32x4::new(-512, 1024, 0, 0).dist_2L1(result) < 8,
                "scale(2,1) * rotate(90) * {:?} => {:?}",
                orig, result);
    }

    #[test]
    fn affine_to_hexagonal() {
        // TODO Shouldn't be necessary; hexagonal_to_cartesian() is
        // unnecessarily complicated
        fn fixup_c(v: i32x4) -> i32x4 {
            i32x4::new(v.extract(0), v.extract(1),
                       -v.extract(0)-v.extract(1), 0)
        }

        let orig = i32x4::new(102400, 25600, 0, 0);
        let rotation = Affine2d::rotate(Wrapping(5461 /* 30 deg */));
        let expected = rotation * orig;
        println!("hexagonal matrix: {:?}", rotation.to_hexagonal());
        let actual = hexgrid::hexagonal_to_cartesian(
            fixup_c(
                rotation.to_hexagonal() *
                hexgrid::cartesian_to_hexagonal(orig)));

        assert!(expected.dist_2L1(actual) < 256,
                "Expected {:?}, got {:?}", expected, actual);
    }

    #[bench]
    fn bench_rotate(b: &mut Bencher) {
        b.iter(|| Affine2d::rotate(black_box(Wrapping(1024))));
    }
}
