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

use super::{Angle,SFrac16,Frac,DEG_90};

const COSINE_TABLE_BITS : usize = 12;
const COSINE_TABLE_SIZE : usize = 1 << COSINE_TABLE_BITS;
const COSINE_TABLE_MASK : usize = COSINE_TABLE_SIZE - 1;
const COSINE_TABLE_SHIFT : usize = 32 - COSINE_TABLE_BITS;

static COSINE_TABLE : [SFrac16; COSINE_TABLE_SIZE] =
    include!("trigtab.rs");

#[inline]
pub fn cos(theta: Angle) -> SFrac16 {
    COSINE_TABLE[(theta >> COSINE_TABLE_SHIFT).0 as usize & COSINE_TABLE_MASK]
}

#[inline]
pub fn sin(theta: Angle) -> SFrac16 {
    cos(theta - DEG_90)
}

#[cfg(test)]
mod tests {
    use std::num::Wrapping;
    use super::super::{SFrac16,Frac,DEG_90,DEG_180,DEG_270};
    use super::*;

    #[test]
    fn simple_cosine() {
        assert_eq!(SFrac16::inverse(1), cos(Wrapping(0)));
        assert_eq!(Frac(0i16), cos(DEG_90));
        assert_eq!(SFrac16::inverse(-1), cos(DEG_180));
        assert_eq!(Frac(0i16), cos(DEG_270));
    }

    #[test]
    fn simple_sin() {
        assert_eq!(Frac(0i16), sin(Wrapping(0)));
        assert_eq!(SFrac16::inverse(1), sin(DEG_90));
        assert_eq!(Frac(0i16), sin(DEG_180));
        assert_eq!(SFrac16::inverse(-1), sin(DEG_270));
    }
}
