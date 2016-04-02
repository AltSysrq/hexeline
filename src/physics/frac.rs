//-
// Copyright (c) 2016, Jason Lingle
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

use std::ops::{Add,Sub,Mul,Neg};
use std::fmt::Debug;

pub trait FractionalBase: Copy + Clone + Debug + PartialEq + Eq + PartialOrd + Ord {
    fn from(x: i64) -> Self;
    fn expand(self) -> i64;
    fn bits() -> i64;
}

impl FractionalBase for i16 {
    fn from(x: i64) -> i16 { x as i16 }
    fn expand(self) -> i64 { self as i64 }
    fn bits() -> i64 { 14 }
}

impl FractionalBase for u32 {
    fn from(x: i64) -> u32 { x as u32 }
    fn expand(self) -> i64 { self as u64 as i64 }
    fn bits() -> i64 { 31 }
}

/// A fixed-point fractional type used to represent values between 0 and 1 or
/// -1 and 1, both inclusive.
///
/// Internally, a fraction is represented by an integer with N bits of
/// precision. If that integer is signed, the fraction represents the range
/// [-1,+1]; if the integer is unsigned, it represents the range [0,+1]. Note
/// that N is generally 1 less than the number of non-sign-bits in the integer
/// so that +1 can be represented exactly. The bits of precision indicates how
/// many steps exist between 0 and +1, regardless of whether negative values
/// are supported.
///
/// When a fraction is mulitplied by an integer, both the underlying
/// representation and the input integer are expanded to 64 bits; those two
/// qwords are multiplied, then the result is bitshifted N bits right and
/// truncated back down to the original integer type.
///
/// Fractions of the same type can be added together, subtracted, and
/// multiplied as if they were regular numeric values. Multiplication does a
/// full fixed-point multiply as described above, such that
/// inverse(2)*inverse(2) == inverse(4).
///
/// Common base types:
///   i16 (SFrac16), 14 bits of precision, -1..+1
///   u32 (UFrac32), 31 bits of precision, 0..+1
#[derive(Copy,Clone,Debug,PartialEq,Eq,PartialOrd,Ord)]
pub struct Frac<T: FractionalBase>(pub T);
pub type SFrac16 = Frac<i16>;
pub type UFrac32 = Frac<u32>;

impl<T: FractionalBase> Frac<T> {
    pub fn inverse(denom: i32) -> Self {
        Self::of(1, denom)
    }

    pub fn of(num: i32, denom: i32) -> Self {
        Frac(T::from(((num as i64) << T::bits()) / (denom as i64)))
    }
}

impl<L,R> Add<Frac<R>> for Frac<L>
where L: FractionalBase + Add<R>,
      R: FractionalBase,
      L::Output: FractionalBase {
    type Output = Frac<L::Output>;

    fn add(self, rhs: Frac<R>) -> Self::Output {
        Frac(self.0 + rhs.0)
    }
}

impl<L,R> Sub<Frac<R>> for Frac<L>
where L: FractionalBase + Sub<R>,
      R: FractionalBase,
      L::Output: FractionalBase {
    type Output = Frac<L::Output>;

    fn sub(self, rhs: Frac<R>) -> Self::Output {
        Frac(self.0 - rhs.0)
    }
}

impl<T: FractionalBase> Mul<Frac<T>> for Frac<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let xlhs = self.0.expand();
        let xrhs = rhs.0.expand();
        let product = xlhs * xrhs;
        Frac(T::from(product >> T::bits()))
    }
}

impl<T> Neg for Frac<T>
where T: FractionalBase + Neg,
      T::Output: FractionalBase {
    type Output = Frac<T::Output>;

    fn neg(self) -> Self::Output {
        Frac(-self.0)
    }
}

impl<T: FractionalBase> Mul<i32> for Frac<T> {
    type Output = i32;

    fn mul(self, rhs: i32) -> i32 {
        let xlhs = self.0.expand();
        let xrhs = rhs as i64;
        let product = xlhs * xrhs;
        (product >> T::bits()) as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn divide_by_2_i16() {
        let frac = SFrac16::inverse(2);
        assert_eq!(2, frac * 4);
        assert_eq!(2, frac * 5);
    }

    #[test]
    fn divide_by_neg_3_i16() {
        let frac = SFrac16::inverse(-3);
        assert_eq!(-2, frac * 6);
        assert_eq!(1, frac * -6);
    }

    #[test]
    fn divide_by_1_i16() {
        let frac = SFrac16::inverse(1);
        assert_eq!(6, frac * 6);
        assert_eq!(-6, frac * -6);
    }

    #[test]
    fn divide_by_neg_1_i16() {
        let frac = SFrac16::inverse(-1);
        assert_eq!(6, frac * -6);
        assert_eq!(-6, frac * 6);
    }

    #[test]
    fn mul_by_two_thirds_u32() {
        let frac = UFrac32::of(2, 3);
        assert_eq!(3, frac * 6);
        assert_eq!(-4, frac * -6);
    }

    #[test]
    fn divide_by_1_u32() {
        let frac = UFrac32::inverse(1);
        assert_eq!(4, frac * 4);
        assert_eq!(-4, frac * -4);
    }

    #[test]
    fn addition() {
        let quarter = UFrac32::inverse(4);
        let half = quarter + quarter;
        assert_eq!(2, half * 4);
    }

    #[test]
    fn subtraction() {
        let half = UFrac32::inverse(2);
        let quarter = UFrac32::inverse(4);
        let q = half - quarter;
        assert_eq!(2, q * 8);
    }

    #[test]
    fn negation() {
        let half = SFrac16::inverse(2);
        let nhalf = -half;
        assert_eq!(-2, nhalf * 4);
    }

    #[test]
    fn multiplication() {
        let half = SFrac16::inverse(2);
        let quarter = half * half;
        assert_eq!(SFrac16::inverse(4), quarter);
        assert_eq!(2, quarter * 8);
    }
}
