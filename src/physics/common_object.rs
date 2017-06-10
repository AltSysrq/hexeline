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
use std::num::Wrapping;

use simd::*;

use physics::Angle;

/// The packed form of the common object data.
///
/// The fields and their representational details are public in order to allow
/// manipulating objects efficiently.
///
/// See `README.org` in this directory for more details on why this is laid out
/// like this.
///
/// Logical field names correspond to those in `UnpackedCommonObject`.
#[derive(Clone, Copy)]
pub struct CommonObject {
    /// Position, collision, and identity data.
    ///
    /// Field 0:
    ///   - [0..31] biased_a
    /// Field 1:
    ///   - [0..31] b
    /// Field 2:
    ///   - [ 0.. 7] rounded_radius
    ///   - [ 8..15] collision_group
    ///   - [16..31] theta
    /// Field 3:
    ///   - [ 0.. 7] object_type
    ///   - [ 8..23] extended_data
    ///   - [24..31] data_dst_size
    ///
    /// Theta is placed at the top of that field so that values can be added to
    /// it without corrupting other data when it wraps around.
    pub p: i32x4,
    /// Dynamics information.
    ///
    /// Field 0:
    ///   - [ 0.. 7] fa
    ///   - [ 8..15] aa
    ///   - [16..31] va
    /// Field 1:
    ///   - [ 0.. 7] fb
    ///   - [ 8..15] ab
    ///   - [16..31] vb
    /// Field 2:
    ///   - [ 0.. 7] ftheta
    ///   - [ 8..15] atheta
    ///   - [16..31] vtheta_x4
    /// Field 3:
    ///   - [ 0.. 7] id_lo
    ///   - [ 8..15] counter_increment
    ///   - [16..23] wakeup_counter
    ///   - [24..31] id_hi
    ///
    /// A and B velocity are aligned so that they can be bit-shifted 16 bits
    /// right and then added to the `a` field. All velocities are in the high
    /// half so that the whole `b` field can be bit-shifted around to
    /// sign-extend acceleration and then added to itself to update the
    /// velocity.
    ///
    /// `counter_increment` and `wakeup_counter` are placed so the counter can
    /// be updated in the same step that updates velocity with acceleration.
    /// Note that this means that the `id` field is awkwardly split in two, and
    /// that `id_hi` gets corrupted when `wakeup_counter` wraps around, which
    /// must be fixed manually.
    pub d: i32x4,
}

/// A representation of `CommonObject` as a normal struct. Use of this is to be
/// avoided except for constructing `CommonObject` and places where performance
/// is unimportant.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnpackedCommonObject {
    /// The A coordinate plus `rounded_radius * ROUNDED_RADIUS_FACTOR`.
    pub biased_a: i32,
    /// The B coordinate.
    pub b: i32,
    /// The rotation.
    pub theta: Angle,
    /// The collision radius of the object (i.e., the maximum distance of any
    /// collideable point from the nominal (A,B) coordinate), divided by
    /// ROUNDED_RADIUS_FACTOR, rounded up.
    pub rounded_radius: u8,
    /// No two objects with the same non-zero collision group will collide with
    /// each other.
    pub collision_group: u8,
    /// The type of this object.
    pub object_type: u8,
    /// If `object_type` identifies an object type with out-of-line data, a
    /// compressed pointer into the extended data heap. Otherwise, some
    /// type-specific data.
    pub extended_data: u16,
    /// If `extended_data` is a pointer to a DST, a type-specific
    /// representation of the length part of the DST pointer. Otherwise, extra
    /// type-specific data.
    pub data_dst_size: u8,
    /// The A velocity.
    pub va: i16,
    /// The A friction, i.e., what fraction of A velocity is preserved every
    /// `FRICTION_TICKS` ticks. 128 means no velocity is lost; 0 means 50% of
    /// velocity is lost.
    pub fa: u8,
    /// The A acceleration (delta `va` per frame).
    pub aa: i8,
    /// The B velocity.
    pub vb: i16,
    /// The B friction, as with `fa`.
    pub fb: u8,
    /// The B acceleration.
    pub ab: i8,
    /// The theta velocity times 4.
    pub vtheta_x4: i16,
    /// The theta friction, as with `fa`.
    pub ftheta: u8,
    /// The theta acceleration, i.e., delta `vtheta_x4` per frame.
    pub atheta: i8,
    /// Counter of time until this object may either invoke extended behaviour
    /// or interact with the static environment. Note that this counts _up_ to
    /// 0 (by wrapping around from 255).
    pub wakeup_counter: u8,
    /// Amount by which `wakeup_counter` is incremented each frame. Currently,
    /// this should always be 1.
    pub wakeup_increment: u8,
    /// The unique id of this object. This is used by graphics and other
    /// processes to track individual objects. A value of 0 indicates an
    /// unidentifiable object. Unidentifiable objects may still have graphics;
    /// they simply cannot have graphical state.
    pub id: u16,
}

macro_rules! common_field {
    ($part:ident:$word:tt[$lo:tt..$hi:tt]: $typ:tt
     $name:ident, $setter:ident, $wither:ident) => {
        common_field!($part:$word[$lo..$hi]: $typ
                      $name, $setter, $wither,
                      $name as $typ,
                      $name as u32 as i32);
    };

    ($part:ident:$word:tt[$lo:tt..$hi:tt]: $typ:tt
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

impl CommonObject {
    common_field!(p:0[ 0..31]: i32 biased_a, set_biased_a, with_biased_a);
    common_field!(p:1[ 0..31]: i32 b, set_b, with_b);
    common_field!(p:2[ 0.. 7]: u8 rounded_radius, set_rounded_radius,
                  with_rounded_radius);
    common_field!(p:2[ 8..15]: u8 collision_group, set_collision_group,
                  with_collision_group);
    common_field!(p:2[16..31]: Angle theta, set_theta, with_theta,
                  Wrapping(theta as i16), theta.0 as i32);
    common_field!(p:3[ 0.. 7]: u8 object_type, set_object_type,
                  with_object_type);
    common_field!(p:3[ 8..23]: u16 extended_data, set_extended_data,
                  with_extended_data);
    common_field!(p:3[24..31]: u8 data_dst_size, set_data_dst_size,
                  with_data_dst_size);

    common_field!(d:0[ 0.. 7]:  u8 fa, set_fa, with_fa);
    common_field!(d:0[ 8..15]:  i8 aa, set_aa, with_aa);
    common_field!(d:0[16..31]: i16 va, set_va, with_va);
    common_field!(d:1[ 0.. 7]:  u8 fb, set_fb, with_fb);
    common_field!(d:1[ 8..15]:  i8 ab, set_ab, with_ab);
    common_field!(d:1[16..31]: i16 vb, set_vb, with_vb);
    common_field!(d:2[ 0.. 7]:  u8 ftheta, set_ftheta, with_ftheta);
    common_field!(d:2[ 8..15]:  i8 atheta, set_atheta, with_atheta);
    common_field!(d:2[16..31]: i16 vtheta_x4, set_vtheta_x4, with_vtheta_x4);
    common_field!(d:3[ 0.. 7]:  u8 id_lo, set_id_lo, with_id_lo);
    common_field!(d:3[ 8..15]:  u8 wakeup_increment, set_wakeup_increment,
                  with_wakeup_increment);
    common_field!(d:3[16..23]:  u8 wakeup_counter, set_wakeup_counter,
                  with_wakeup_counter);
    common_field!(d:3[24..31]:  u8 id_hi, set_id_hi, with_id_hi);

    pub fn unpack(self) -> UnpackedCommonObject {
        UnpackedCommonObject {
            biased_a: self.biased_a(),
            b: self.b(),
            rounded_radius: self.rounded_radius(),
            collision_group: self.collision_group(),
            theta: self.theta(),
            object_type: self.object_type(),
            extended_data: self.extended_data(),
            data_dst_size: self.data_dst_size(),
            fa: self.fa(),
            aa: self.aa(),
            va: self.va(),
            fb: self.fb(),
            ab: self.ab(),
            vb: self.vb(),
            vtheta_x4: self.vtheta_x4(),
            atheta: self.atheta(),
            ftheta: self.ftheta(),
            wakeup_counter: self.wakeup_counter(),
            wakeup_increment: self.wakeup_increment(),
            id: (self.id_lo() as u16) | ((self.id_hi() as u16) << 8),
        }
    }
}

impl UnpackedCommonObject {
    #[inline]
    pub fn pack(self) -> CommonObject {
        CommonObject {
            p: i32x4::splat(0),
            d: i32x4::splat(0),
        }
            .with_biased_a(self.biased_a)
            .with_b(self.b)
            .with_rounded_radius(self.rounded_radius)
            .with_collision_group(self.collision_group)
            .with_theta(self.theta)
            .with_object_type(self.object_type)
            .with_extended_data(self.extended_data)
            .with_data_dst_size(self.data_dst_size)
            .with_fa(self.fa)
            .with_aa(self.aa)
            .with_va(self.va)
            .with_fb(self.fb)
            .with_ab(self.ab)
            .with_vb(self.vb)
            .with_ftheta(self.ftheta)
            .with_atheta(self.atheta)
            .with_vtheta_x4(self.vtheta_x4)
            .with_wakeup_counter(self.wakeup_counter)
            .with_wakeup_increment(self.wakeup_increment)
            .with_id_lo(self.id as u8)
            .with_id_hi((self.id >> 8) as u8)
    }
}

impl fmt::Debug for CommonObject {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(&self.unpack(), f)
    }
}

#[cfg(test)]
mod test {
    use std::{i8, u8, i16, u16, i32, u32};
    use std::num::Wrapping;

    use super::*;

    #[test]
    fn packed_object_stores_all_fields_correctly() {
        const I32_BOUNDARIES: &[i32] = &[0, i32::MIN, i32::MAX];
        const U32_BOUNDARIES: &[u32] = &[0, u32::MAX];
        const I16_BOUNDARIES: &[i16] = &[0, i16::MIN, i16::MAX];
        const U16_BOUNDARIES: &[u16] = &[0, u16::MAX];
        const I8_BOUNDARIES:  &[i8 ] = &[0, i8::MIN, i8::MAX];
        const U8_BOUNDARIES:  &[u8 ] = &[0, u8::MAX];
        const ID_BOUNDARIES:  &[u16] = &[0, 128, 255, 256, u16::MAX];

        // Brute-force test of all boundary values
        for &biased_a in I32_BOUNDARIES {
        for &b in I32_BOUNDARIES {
        for theta in I16_BOUNDARIES.iter().map(|&t| Wrapping(t)) {
        for &rounded_radius in U8_BOUNDARIES {
        for &collision_group in U8_BOUNDARIES {
        for &object_type in U8_BOUNDARIES {
        for &extended_data in U16_BOUNDARIES {
        for &data_dst_size in U8_BOUNDARIES {
        for &va in I16_BOUNDARIES {
        for &fa in U8_BOUNDARIES {
        for &aa in I8_BOUNDARIES {
        for &vb in I16_BOUNDARIES {
        for &fb in U8_BOUNDARIES {
        for &ab in I8_BOUNDARIES {
        for &vtheta_x4 in I16_BOUNDARIES {
        for &ftheta in U8_BOUNDARIES {
        for &atheta in I8_BOUNDARIES {
        for &wakeup_counter in U8_BOUNDARIES {
        for &wakeup_increment in U8_BOUNDARIES {
        for &id in ID_BOUNDARIES {
            let orig = UnpackedCommonObject {
                biased_a, b, theta, rounded_radius,
                collision_group, object_type, extended_data,
                data_dst_size,
                va, fa, aa, vb, fb, ab, vtheta_x4, ftheta, atheta,
                wakeup_counter, wakeup_increment, id,
            };
            let packed = orig.pack();
            let unpacked = packed.unpack();
            assert_eq!(orig, unpacked);
        }}}}}}}}}}}}}}}}}}}}
    }
}
