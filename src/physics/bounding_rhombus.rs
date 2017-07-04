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

use std::fmt;

use simd::*;
use simdext::*;

use physics::coords::*;

/// Represents an axis-aligned region of hexagonal space.
///
/// This is the hexagonal equivalent of a bounding rectangle.
///
/// The bounding area is inclusive on both bounds.
///
/// Internally, this is stored as the *negated* minimum (A,B) coordinates in
/// the first half of a SIMD word and with the (non-negated) maximum (A,B)
/// coordinates in the second half. This allows for very efficient unions and
/// overlap checks.
#[derive(Clone, Copy)]
pub struct BoundingRhombus(i32x4);

/// A `BoundingRhombus` which has been "inverted" for fast overlap checks.
#[derive(Clone, Copy)]
pub struct InvertedBoundingRhombus(i32x4);

impl fmt::Debug for BoundingRhombus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let dim = self.dim();
        write!(f, "BoundingRhombus({:?}..{:?})",
               dim.fst_dual(), dim.snd_dual())
    }
}

impl BoundingRhombus {
    /// Create a `BoundingRhombus` centred on `point` and extending `radius`
    /// distance along the A and B axes.
    ///
    /// Beware that `radius` does not correspond to the radius of the bounding
    /// circle, since the unit circle has an L2 distance from the centre
    /// greater than 1 in hexagonal space.
    #[inline(always)]
    pub fn around(point: Vhs, radius: i32) -> Self {
        let point = point.repr();
        let neg = i32x4::splat(0) - point;
        let centre = i32x4::new(neg.extract(0), neg.extract(1),
                                point.extract(0), point.extract(1));
        BoundingRhombus(centre + i32x4::splat(radius))
    }

    /// Create a `BoundingRhombus` covering the area between the two points
    /// given.
    ///
    /// Point 1 must have coordinates less than or equal to those of point 2.
    #[inline(always)]
    pub fn of(area: Vhl) -> Self {
        debug_assert!(area.a1() <= area.a2() &&
                      area.b1() <= area.b2(),
                      "Improperly ordered rhombus area: {:?}", area);
        let area = area.repr();
        let neg = i32x4::splat(0) - area;
        BoundingRhombus(i32x4::new(neg.extract(0), neg.extract(1),
                                   area.extract(2), area.extract(3)))
    }

    /// Return the two boundary points of this `BoundingRhombus`, low point
    /// first.
    #[inline(always)]
    pub fn dim(self) -> Vhl {
        let neg = i32x4::splat(0) - self.0;
        Vhl(neg.extract(0), neg.extract(1),
            self.0.extract(2), self.0.extract(3))
    }

    /// Create a `BoundingRhombus` which covers at least all space that is
    /// covered by either `self` or `other`.
    #[inline(always)]
    pub fn union(self, other: BoundingRhombus) -> Self {
        BoundingRhombus(self.0.nsw_max(other.0))
    }

    /// Preprocess this `BoundingRhombus` to make it ready for `overlaps()`
    /// tests.
    #[inline(always)]
    pub fn invert(self) -> InvertedBoundingRhombus {
        // See `overlaps` for what this means
        InvertedBoundingRhombus((i32x4::splat(0) - self.0).shuf(2, 3, 0, 1))
    }

    /// Return whether `self` and `other` overlap at any point.
    #[inline(always)]
    pub fn overlaps(self, other: InvertedBoundingRhombus) -> bool {
        // Two intervals (a0..a1),(b0..b1) overlap iff
        //   a0 <= b1 && b0 <= a1
        // We store the lower bounds negated, which makes it initially more
        // difficult:
        //   -(-a0) <= b1 && -(-b0) <= a1
        // However, we can simply negate both B terms and thus handle the
        // negation with one operation:
        //   (-a0) >= -b1 && -(-b0) <= a1
        //
        // So this means we have two comparisons per axis, or four in total.
        // This means we can do all four at once provided we can use the same
        // operation for all of them. And in fact we can, by simply ensuring
        // the b term is on the right-hand side:
        //   (-a0) >= -b1 && a1 >= -(-b0)

        let diff = self.0 - other.0;
        // If a >= b, that lane will have the sign bit clear. So all four sign
        // bits must be clear for there to be overlap.
        !diff.any_sign_bit()
    }
}

#[cfg(test)]
mod test {
    use test::{Bencher, black_box};

    use super::*;

    #[bench]
    fn bench_around(b: &mut Bencher) {
        b.iter(|| BoundingRhombus::around(
            black_box(Vhs(3, 4)), black_box(5)))
    }

    #[bench]
    fn bench_union(b: &mut Bencher) {
        b.iter(|| black_box(
            BoundingRhombus::around(Vhs(3, 4), 5)).union(
            black_box(
                BoundingRhombus::around(Vhs(5, 6), 7))))
    }

    #[bench]
    fn bench_overlaps(b: &mut Bencher) {
        b.iter(|| black_box(
            BoundingRhombus::around(Vhs(3, 4), 5)).overlaps(
            black_box(
                BoundingRhombus::around(Vhs(5, 6), 7)).invert()))
    }

    #[test]
    fn test_overlap() {
        struct Case(i32, i32, i32, i32,
                    i32, i32, i32, i32,
                    bool);

        let cases = [
            // Identical rhombi
            Case(-10, -10, 10, 10,
                 -10, -10, 10, 10,
                 true),
            // No axis overlap, various relative positions
            Case(-10, -10, 10, 10,
                 -20, -20, -15, -15,
                 false),
            Case(-10, -10, 10, 10,
                 -20, 15, -15, 20,
                 false),
            Case(-10, -10, 10, 10,
                 15, -20, 20, -15,
                 false),
            Case(-10, -10, 10, 10,
                 15, 15, 20, 20,
                 false),
            // Overlap on one axis, none on the other
            Case(-10, -10, 10, 10,
                 -15, -20, -5, -15,
                 false),
            Case(-10, -10, 10, 10,
                 -5, -20, 5, -15,
                 false),
            Case(-10, -10, 10, 10,
                 5, -20, 15, -15,
                 false),
            Case(-10, -10, 10, 10,
                 -20, -20, -15, -5,
                 false),
            Case(-10, -10, 10, 10,
                 -20, -5, -15, 5,
                 false),
            Case(-10, -10, 10, 10,
                 -20, 5, -15, 15,
                 false),
            // Corner contained in other
            Case(-10, -10, 10, 10,
                 -20, -20, 0, 0,
                 true),
            Case(-10, -10, 10, 10,
                 -20, 0, 0, 20,
                 true),
            Case(-10, -10, 10, 10,
                 0, -20, 20, 0,
                 true),
            Case(-10, -10, 10, 10,
                 0, 0, 20, 20,
                 true),
            // Two corners contained in other
            Case(-10, -10, 10, 10,
                 -20, -5, 0, 5,
                 true),
            Case(-10, -10, 10, 10,
                 -5, -20, 5, 0,
                 true),
            Case(-10, -10, 10, 10,
                 0, -5, 20, 5,
                 true),
            Case(-10, -10, 10, 10,
                 -5, 0, 5, 20,
                 true),
            // One box completely contained in other
            Case(-10, -10, 10, 10,
                 -20, -20, 20, 20,
                 true),
            Case(-10, -10, 10, 10,
                 -5, -5, 5, 5,
                 true),
            // Edges only overlapping
            Case(-10, -10, 10, 10,
                 -20, -20, 20, -10,
                 true),
            Case(-10, -10, 10, 10,
                 -20, 10, 20, 20,
                 true),
            Case(-10, -10, 10, 10,
                 -20, -20, -10, 20,
                 true),
            Case(-10, -10, 10, 10,
                 10, -20, 20, 20,
                 true),
            // Only corner touching
            Case(-10, -10, 10, 10,
                 -20, -20, -10, -10, // Touch at (-10,-10)
                 true),
            Case(-10, -10, 10, 10,
                 -10, 10, 10, 20, // Touch at (-10, 10)
                 true),
            Case(-10, -10, 10, 10,
                 10, -10, 20, 10, // Touch at (10, -10)
                 true),
            Case(-10, -10, 10, 10,
                 10, 10, 20, 20, // Touch at (10,10)
                 true),
        ];

        for case in &cases {
            let a = BoundingRhombus::of(Vhl(case.0, case.1, case.2, case.3));
            let b = BoundingRhombus::of(Vhl(case.4, case.5, case.6, case.7));
            assert_eq!(case.8, a.overlaps(b.invert()),
                       "Incorrect result for overlapping of {:?} and {:?}",
                       a, b);
            assert_eq!(case.8, b.overlaps(a.invert()),
                       "Non-transitive result for overlapping \
                        {:?} and {:?}", a, b);
        }
    }

    #[test]
    fn test_union() {
        let a = BoundingRhombus::of(Vhl(-20, -20, 10, 10));
        let b = BoundingRhombus::of(Vhl(-10, -10, 20, 20));

        let u = a.union(b).dim();
        assert_eq!(-20, u.a1());
        assert_eq!(-20, u.b1());
        assert_eq!(20, u.a2());
        assert_eq!(20, u.b2());
        let u = b.union(a).dim();
        assert_eq!(-20, u.a1());
        assert_eq!(-20, u.b1());
        assert_eq!(20, u.a2());
        assert_eq!(20, u.b2());

        let a = BoundingRhombus::of(Vhl(-20, -10, 10, 20));
        let b = BoundingRhombus::of(Vhl(-10, -20, 20, 10));

        let u = a.union(b).dim();
        assert_eq!(-20, u.a1());
        assert_eq!(-20, u.b1());
        assert_eq!(20, u.a2());
        assert_eq!(20, u.b2());
        let u = b.union(a).dim();
        assert_eq!(-20, u.a1());
        assert_eq!(-20, u.b1());
        assert_eq!(20, u.a2());
        assert_eq!(20, u.b2());
    }
}
