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
use std::i32;

use simd::*;
use simdext::*;

use physics::Angle;
use physics::coords::Vhs;

pub const ROUNDED_RADIUS_SHIFT: u32 = 8;
pub const ROUNDED_RADIUS_FACTOR: u32 = 1 << ROUNDED_RADIUS_SHIFT;

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
    ///   - [ 0..31] a
    /// Field 1:
    ///   - [ 0..31] b
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
    ///   - [ 8..15] aax16
    ///   - [16..31] vax4
    /// Field 1:
    ///   - [ 0.. 7] fb
    ///   - [ 8..15] abx16
    ///   - [16..31] vbx4
    /// Field 2:
    ///   - [ 0.. 7] ftheta
    ///   - [ 8..15] athetax16
    ///   - [16..31] vthetax4
    /// Field 3:
    ///   - [ 0.. 7] id_lo
    ///   - [ 8..15] counter_increment
    ///   - [16..23] wakeup_counter
    ///   - [24..31] id_hi
    ///
    /// A and B velocity are aligned so that they can be bit-shifted 16 bits
    /// right and then added to the position field. All velocities are in the
    /// high half so that the whole `b` field can be bit-shifted around to
    /// sign-extend acceleration and then added to itself to update the
    /// velocity.
    ///
    /// Nominal velocities and accelerations are multiplied by 4 and 16,
    /// respectively, to give better precision and a more useful range:
    ///
    /// | Coordinate            | Internal Unit         | Natural Unit          |
    /// |=======================|=======================|=======================|
    /// | Min A velocity / 1    | 1 space/tick          | 3 px/sec              |
    /// | Max A velocity / 1    | 32768 space/tick      | 50 screens/sec        |
    /// | Min A accel / 1       | 1 space/tick/tick     | 0.2 screens/sec/sec   |
    /// | Max A accel / 1       | 128 space/tick/tick   | 19 screens/sec/sec    |
    /// | Min theta velocity/1  | 1 rot/tick            | 0.5 deg/sec           |
    /// | Max theta velocity/1  | 32768 rot/tick        | 50 spin/sec           |
    /// | Min theta accel/1     | 1 rot/tick/tick       | 20 deg/sec/sec        |
    /// | Max theta accel/1     | 128 rot/tick/tick     | 19.5 spin/sec/sec     |
    /// |                       |                       |                       |
    /// | Min A velocity / 4    | 0.25 space/tick       | 0.73 px/sec           |
    /// | Max A velocity / 4    | 8192 space/tick       | 25 screens/sec        |
    /// | Min A accel / 16      | 1/16 space/tick/tick  | 18 px/sec/sec         |
    /// | Max A accel / 16      | 8 space/tick/tick     | 1.2 screens/sec/sec   |
    /// | Min theta velocity/4  | 0.25 rot/tick         | 0.137 deg/sec         |
    /// | Max theta velocity/4  | 8192 rot/tick         | 12.5 spins/sec        |
    /// | Min theta accel / 16  | 1/16 rot/tick/tick    | 3.4 deg/sec/sec       |
    /// | Max theta accel / 16  | 8 rot/tick/tick       | 439 deg/sec/sec       |
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct UnpackedCommonObject {
    /// The A coordinate.
    pub a: i32,
    /// The B coordinate.
    pub b: i32,
    /// The rotation.
    pub theta: Angle,
    /// The collision radius of the object (i.e., the maximum displacement
    /// along any *one* hexagonal axis to any collideable part of the object),
    /// divided by ROUNDED_RADIUS_FACTOR, rounded up.
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
    /// The A velocity multiplied by 4.
    pub vax4: i16,
    /// The A friction, i.e., what fraction of A velocity is preserved every
    /// tick. This is implicitly biased by 65281, with a final value of 65536
    /// indicating no velocity is lost, and 0 indicating retaining around 68%
    /// of velocity per second.
    pub fa: u8,
    /// The A acceleration (delta `va` per frame) multiplied by 16.
    pub aax16: i8,
    /// The B velocity multiplied by 4.
    pub vbx4: i16,
    /// The B friction, as with `fa`.
    pub fb: u8,
    /// The B acceleration multiplied by 16.
    pub abx16: i8,
    /// The theta velocity times 4.
    pub vthetax4: i16,
    /// The theta friction, as with `fa`.
    pub ftheta: u8,
    /// The theta acceleration multiplied by 16.
    pub athetax16: i8,
    /// Counter of time until this object may either invoke extended behaviour
    /// or interact with the static environment. Note that this counts _up_ to
    /// 255.
    pub wakeup_counter: u8,
    /// Amount by which `wakeup_counter` is incremented each frame, divided by
    /// 4. Currently, this should always be 4.
    pub wakeup_increment: u8,
    /// The unique id of this object. This is used by graphics and other
    /// processes to track individual objects. A value of 0 indicates an
    /// unidentifiable object. Unidentifiable objects may still have graphics;
    /// they simply cannot have graphical state.
    pub id: u16,
}

impl CommonObject {
    packed_field!(p:0[ 0..31]: i32 a, set_a, with_a);
    packed_field!(p:1[ 0..31]: i32 b, set_b, with_b);
    packed_field!(p:2[ 0.. 7]: u8 rounded_radius, set_rounded_radius,
                  with_rounded_radius);
    packed_field!(p:2[ 8..15]: u8 collision_group, set_collision_group,
                  with_collision_group);
    packed_field!(p:2[16..31]: Angle theta, set_theta, with_theta,
                  Wrapping(theta as i16), theta.0 as i32);
    packed_field!(p:3[ 0.. 7]: u8 object_type, set_object_type,
                  with_object_type);
    packed_field!(p:3[ 8..23]: u16 extended_data, set_extended_data,
                  with_extended_data);
    packed_field!(p:3[24..31]: u8 data_dst_size, set_data_dst_size,
                  with_data_dst_size);

    packed_field!(d:0[ 0.. 7]:  u8 fa   , set_fa   , with_fa   );
    packed_field!(d:0[ 8..15]:  i8 aax16, set_aax16, with_aax16);
    packed_field!(d:0[16..31]: i16 vax4 , set_vax4 , with_vax4 );
    packed_field!(d:1[ 0.. 7]:  u8 fb   , set_fb   , with_fb   );
    packed_field!(d:1[ 8..15]:  i8 abx16, set_abx16, with_abx16);
    packed_field!(d:1[16..31]: i16 vbx4 , set_vbx4 , with_vbx4 );
    packed_field!(d:2[ 0.. 7]:  u8 ftheta   , set_ftheta   , with_ftheta   );
    packed_field!(d:2[ 8..15]:  i8 athetax16, set_athetax16, with_athetax16);
    packed_field!(d:2[16..31]: i16 vthetax4 , set_vthetax4 , with_vthetax4 );
    packed_field!(d:3[ 0.. 7]:  u8 id_lo, set_id_lo, with_id_lo);
    packed_field!(d:3[ 8..15]:  u8 wakeup_increment, set_wakeup_increment,
                  with_wakeup_increment);
    packed_field!(d:3[16..23]:  u8 wakeup_counter, set_wakeup_counter,
                  with_wakeup_counter);
    packed_field!(d:3[24..31]:  u8 id_hi, set_id_hi, with_id_hi);

    pub fn id(self) -> u16 {
        (self.id_lo() as u16) | ((self.id_hi() as u16) << 8)
    }

    pub fn unpack(self) -> UnpackedCommonObject {
        UnpackedCommonObject {
            a: self.a(),
            b: self.b(),
            rounded_radius: self.rounded_radius(),
            collision_group: self.collision_group(),
            theta: self.theta(),
            object_type: self.object_type(),
            extended_data: self.extended_data(),
            data_dst_size: self.data_dst_size(),
            fa: self.fa(),
            aax16: self.aax16(),
            vax4: self.vax4(),
            fb: self.fb(),
            abx16: self.abx16(),
            vbx4: self.vbx4(),
            vthetax4: self.vthetax4(),
            athetax16: self.athetax16(),
            ftheta: self.ftheta(),
            wakeup_counter: self.wakeup_counter(),
            wakeup_increment: self.wakeup_increment(),
            id: self.id(),
        }
    }

    /// Return the position of this object as a vector.
    pub fn pos(&self) -> Vhs {
        Vhs::from_repr(self.p)
    }

    /// Advance this object forward one tick, considering only the common data.
    ///
    /// `tick_mod_4` is the lower 2 bits of the current tick number, used to
    /// modulate the `x4` and `x16` fields.
    pub fn tick(self, tick_mod_4: i32) -> Self {
        let velocity: i32x4 = self.d >> 16;
        let acceleration: i32x4 = self.d << 16 >> 24;

        let mod_velocity: i32x4 = (velocity + i32x4::splat(tick_mod_4)) >> 2;
        let real_velocity = i32x4::new(
            mod_velocity.extract(0), mod_velocity.extract(1),
            mod_velocity.extract(2) << 16, 0);

        let mod_acceleration: i32x4 =
            (acceleration + i32x4::splat(tick_mod_4)) >> 2;

        let new_pos = self.p + real_velocity;
        // This also increments wakeup_counter
        let new_velocity = velocity + mod_acceleration;
        // Clamp velocity to legal values.
        //
        // This needs to happen *before* the friction multiplication to ensure
        // it does not overflow.
        let new_velocity = new_velocity
            .nsw_clamp3(i32x4::new(-32768, -32768, -32768, i32::MIN),
                        i32x4::new(32767, 32767, 32767, i32::MAX));
        let friction = (self.d & i32x4::new(255, 255, 255, 0)) +
            i32x4::new(65281, 65281, 65281, 65536);
        let new_velocity: i32x4 = new_velocity * friction >> 16;

        let new_dynamics =
            (new_velocity << 16) |
            (self.d & i32x4::splat(0xFFFF));

        CommonObject { p: new_pos, d: new_dynamics }
    }

    /// Translate the given l2 distance into a (non-rounded) object radius.
    #[inline(always)]
    pub fn radius_of_l2(l2: u32) -> u32 {
        // The maximum contribution of any single hexagonal axis to the L2
        // distance is given by
        //
        //       |a|
        // --------------
        // sqrt(a²+b²+c²)
        //
        // Since c=-a-b,
        //
        //       |a|
        // -----------------
        // sqrt(2a²+2b²-2ab)
        //
        // The maximum value of this function is sqrt(2/3)
        ((l2 as u64 * 53570) >> 16) as u32
    }

    /// Format `radius` as a suitable value for `rounded_radius`.
    #[inline(always)]
    pub fn round_radius(radius: u32) -> u8 {
        debug_assert!(radius <= 255 << ROUNDED_RADIUS_SHIFT);
        ((radius + ROUNDED_RADIUS_FACTOR - 1) >> ROUNDED_RADIUS_SHIFT) as u8
    }
}

impl UnpackedCommonObject {
    #[inline]
    pub fn pack(self) -> CommonObject {
        CommonObject {
            p: i32x4::splat(0),
            d: i32x4::splat(0),
        }
            .with_a(self.a)
            .with_b(self.b)
            .with_rounded_radius(self.rounded_radius)
            .with_collision_group(self.collision_group)
            .with_theta(self.theta)
            .with_object_type(self.object_type)
            .with_extended_data(self.extended_data)
            .with_data_dst_size(self.data_dst_size)
            .with_fa(self.fa)
            .with_aax16(self.aax16)
            .with_vax4(self.vax4)
            .with_fb(self.fb)
            .with_abx16(self.abx16)
            .with_vbx4(self.vbx4)
            .with_ftheta(self.ftheta)
            .with_athetax16(self.athetax16)
            .with_vthetax4(self.vthetax4)
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
    use std::cmp::{max, min};
    use std::{i8, u8, i16, u16, i32, u32};
    use std::num::Wrapping;
    use test::{black_box, Bencher};

    use proptest;
    use proptest::strategy::Strategy;

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
        for &a in I32_BOUNDARIES {
        for &b in I32_BOUNDARIES {
        for theta in I16_BOUNDARIES.iter().map(|&t| Wrapping(t)) {
        for &rounded_radius in U8_BOUNDARIES {
        for &collision_group in U8_BOUNDARIES {
        for &object_type in U8_BOUNDARIES {
        for &extended_data in U16_BOUNDARIES {
        for &data_dst_size in U8_BOUNDARIES {
        for &vax4 in I16_BOUNDARIES {
        for &fa in U8_BOUNDARIES {
        for &aax16 in I8_BOUNDARIES {
        for &vbx4 in I16_BOUNDARIES {
        for &fb in U8_BOUNDARIES {
        for &abx16 in I8_BOUNDARIES {
        for &vthetax4 in I16_BOUNDARIES {
        for &ftheta in U8_BOUNDARIES {
        for &athetax16 in I8_BOUNDARIES {
        for &wakeup_counter in U8_BOUNDARIES {
        for &wakeup_increment in U8_BOUNDARIES {
        for &id in ID_BOUNDARIES {
            let orig = UnpackedCommonObject {
                a, b, theta, rounded_radius,
                collision_group, object_type, extended_data,
                data_dst_size,
                vax4, fa, aax16, vbx4, fb, abx16, vthetax4, ftheta, athetax16,
                wakeup_counter, wakeup_increment, id,
            };
            let packed = orig.pack();
            let unpacked = packed.unpack();
            assert_eq!(orig, unpacked);
        }}}}}}}}}}}}}}}}}}}}
    }

    #[test]
    fn no_vtheta_wraparound() {
        let start = UnpackedCommonObject {
            vthetax4: -32740,
            ftheta: 252,
            athetax16: -120,
            .. UnpackedCommonObject::default()
        }.pack();
        let result = start.tick(0);
        assert!(result.vthetax4() < 0,
                "vthetax4 wrapped around: {} => {}",
                start.vthetax4(), result.vthetax4());
    }

    #[test]
    fn no_vtheta_wraparound_2() {
        let start = UnpackedCommonObject {
            vthetax4: 32738,
            athetax16: 118,
            ftheta: 255,
            wakeup_increment: 4,
            .. UnpackedCommonObject::default()
        }.pack();
        let result = start.tick(2);
        assert_eq!(32767, result.vthetax4());
    }

    #[test]
    fn id_not_clobbered_by_tick() {
        let start = UnpackedCommonObject {
            wakeup_increment: 4,
            .. UnpackedCommonObject::default()
        }.pack();
        let result = start.tick(0);
        assert_eq!(0, result.id());
    }

    proptest! {
        #![proptest_config(proptest::test_runner::Config {
            cases: 2048,
            .. proptest::test_runner::Config::default()
        })]

        #[test]
        fn common_object_tick(
            a in -0x10000000i32..0x10000000,
            b in -0x10000000i32..0x10000000,
            theta in proptest::num::i16::ANY.prop_map(Wrapping),
            rounded_radius in proptest::num::u8::ANY,
            collision_group in proptest::num::u8::ANY,
            object_type in proptest::num::u8::ANY,
            extended_data in proptest::num::u16::ANY,
            data_dst_size in proptest::num::u8::ANY,
            vax4 in proptest::num::i16::ANY,
            fa in proptest::num::u8::ANY,
            aax16 in proptest::num::i8::ANY,
            vbx4 in proptest::num::i16::ANY,
            fb in proptest::num::u8::ANY,
            abx16 in proptest::num::i8::ANY,
            vthetax4 in proptest::num::i16::ANY,
            ftheta in proptest::num::u8::ANY,
            athetax16 in proptest::num::i8::ANY,
            wakeup_counter in 0u8..255u8, // Can't wrap around
            id in proptest::num::u16::ANY,
            tick_mod_4 in 0i32..4i32
        ) {
            let start = UnpackedCommonObject {
                a, b, theta, rounded_radius, collision_group,
                object_type, extended_data, data_dst_size,
                vax4, fa, aax16, vbx4, fb, abx16, vthetax4, ftheta,
                athetax16, wakeup_counter, id,
                wakeup_increment: 4,
            }.pack();
            let result = start.tick(tick_mod_4);

            // Constant properties don't change
            assert_eq!(rounded_radius, result.rounded_radius());
            assert_eq!(collision_group, result.collision_group());
            assert_eq!(object_type, result.object_type());
            assert_eq!(extended_data, result.extended_data());
            assert_eq!(data_dst_size, result.data_dst_size());
            assert_eq!(fa, result.fa());
            assert_eq!(aax16, result.aax16());
            assert_eq!(fb, result.fb());
            assert_eq!(abx16, result.abx16());
            assert_eq!(ftheta, result.ftheta());
            assert_eq!(athetax16, result.athetax16());
            assert_eq!(id, result.id());
            assert_eq!(4, result.wakeup_increment());

            // Position adjusted by velocity or velocity + 1
            // Note that we need to use >>2 and not /4 because they give
            // different results for negative values (-15841>>2 = -3961,
            // -15841/4 = -3960).
            assert!(a + (vax4 as i32 >> 2) == result.a() ||
                    a + (vax4 as i32 >> 2) + 1 == result.a(),
                    "a ({}) + vax4 ({}) -> {}",
                    a, vax4, result.a());
            assert!(b + (vbx4 as i32 >> 2) == result.b() ||
                    b + (vbx4 as i32 >> 2) + 1 == result.b(),
                    "b ({}) + vbx4 ({}) -> {}",
                    b, vbx4, result.b());
            assert!(theta + Wrapping(vthetax4 as i16 >> 2) == result.theta() ||
                    theta + Wrapping(vthetax4 as i16 >> 2) +
                    Wrapping(1) == result.theta(),
                    "theta ({}) + vthetax4 ({}) -> {}",
                    theta, vthetax4, result.theta());

            // The interaction between acceleration and friction on velocity is
            // hard to test at the same time; they are tested in separate
            // tests. But do test that we didn't see any wrap-around of any of
            // the velocities.
            assert!((vax4 as i32 - result.vax4() as i32).abs() < 16384,
                    "Velocity wrap-around for A: {} => {}",
                    vax4, result.vax4());
            assert!((vbx4 as i32 - result.vbx4() as i32).abs() < 16384,
                    "Velocity wrap-around for B: {} => {}",
                    vbx4, result.vbx4());
            assert!((vthetax4 as i32 - result.vthetax4() as i32).abs() < 16384,
                    "Velocity wrap-around for theta: {} => {}",
                    vthetax4, result.vthetax4());

            assert_eq!(wakeup_counter + 1, result.wakeup_counter());

            // Test friction
            let no_accel_start = start
                .with_aax16(0)
                .with_abx16(0)
                .with_athetax16(0);
            let no_accel_result = no_accel_start.tick(tick_mod_4);
            macro_rules! check_friction {
                ($v:ident, $f:ident) => {
                    if 255 == $f {
                        assert_eq!($v, no_accel_result.$v());
                    } else {
                        assert!((no_accel_result.$v() as i32).abs() <=
                                ($v as i32).abs() &&
                                (no_accel_result.$v() as i32).abs() >=
                                ($v as i32).abs()/2,
                                "{} ({}) * friction ({}) -> {}",
                                stringify!($v), $v, $f, no_accel_result.$v());
                    }
                }
            };
            check_friction!(vax4, fa);
            check_friction!(vbx4, fb);
            check_friction!(vthetax4, ftheta);

            // Test acceleration
            let no_friction_start = start
                .with_fa(255)
                .with_fb(255)
                .with_ftheta(255);
            let no_friction_result = no_friction_start.tick(tick_mod_4);
            macro_rules! check_accel {
                ($v:ident, $a:ident) => {
                    let new_v = max(-32768, min(32767, $v as i32 + ($a as i32 >>2)));
                    assert!(no_friction_result.$v() as i32 == new_v ||
                            no_friction_result.$v() as i32 == new_v + 1,
                            "{} ({}) + {} -> {}",
                            stringify!($v), $v, $a, no_friction_result.$v());
                }
            };
            check_accel!(vax4, aax16);
            check_accel!(vbx4, abx16);
            check_accel!(vthetax4, athetax16);
        }
    }

    #[bench]
    fn bench_common_object_tick(b: &mut Bencher) {
        b.iter(|| black_box(CommonObject {
            p: i32x4::splat(0),
            d: i32x4::splat(0),
        }).tick(black_box(1)));
    }
}
