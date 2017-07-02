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

#[macro_export]
macro_rules! packed_field {
    ($part:tt:$word:tt[$lo:tt..$hi:tt]: $typ:tt
     $name:ident, $setter:ident, $wither:ident) => {
        packed_field!($part:$word[$lo..$hi]: $typ
                      $name, $setter, $wither,
                      $name as $typ,
                      $name as u32 as i32);
    };

    ($part:tt:$word:tt[$lo:tt..$hi:tt]: $typ:tt
     $name:ident, $setter:ident, $wither:ident,
     $from_i32:expr, $to_i32:expr) => {
        #[inline]
        pub fn $name(self) -> $typ {
            let $name = self.$part.extract($word) >> $lo;
            $from_i32
        }

        #[inline]
        pub fn $setter(&mut self, $name: $typ) {
            let v: i32 = $to_i32;
            let mask = (((1i64 << ($hi+1)) - 1) ^ ((1i64 << $lo) - 1)) as i32;
            let mut word = self.$part.extract($word);
            word &= !mask;
            word |= (v << $lo) & mask;
            self.$part = self.$part.replace($word, word);
        }

        #[inline]
        pub fn $wither(mut self, $name: $typ) -> Self {
            self.$setter($name);
            self
        }
    };
}

