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

//! `ExtendedObject` types used only for testing.

use std::num::Wrapping;

use simd::i32x4;
use smallvec::SmallVec;

use super::{ExtendedObject, TYPE_DEBUG_PLAYER, TYPE_DEBUG_PROJECTILE};
use physics::common_object::*;
use physics::composite::*;
use physics::coords::*;
use physics::event::*;
use physics::snapshot::*;
use physics::xform::Affine2dH;

/// A point-particle which ceases to exist the instant it undergoes a
/// non-trivial update.
eo_deftype!(DebugProjectile {
    fn update(&self, _pipeline: &mut SnapshotUpdatePipeline,
              _spawn: &mut SpawnList,
              _collisions: &[ImpendingCollision],
              _snapshot: &Snapshot) {
        // Do nothing. Projectile ceases to exist.
    }
});

/// A composite which tracks input events and responds to some of them.
eo_deftype!(DebugPlayer {
    fn update(&self, pipeline: &mut SnapshotUpdatePipeline,
              spawn: &mut SpawnList,
              _collisions: &[ImpendingCollision],
              snapshot: &Snapshot) {
        let mut self_obj = self.0;
        // Reset the wakeup counter so we can get events next frame
        self_obj.set_wakeup_counter(255);

        let self_comp = unsafe {
            CompositeObject::wrap(snapshot.decompress_obj(self_obj))
        };
        let old_ctrl = ConstantControlFlags::from_bits_truncate(
            self_comp.tail()[0].extract(0) as u32);
        let old_rotate = self_comp.tail()[0].extract(1) as i16;
        let last_fire = self_comp.tail()[0].extract(2) as u32;

        let mut new_ctrl = old_ctrl;
        let mut new_rotate = old_rotate;
        let mut new_fire = last_fire;

        for event in pipeline.poll_events(self_obj.id()) {
            match *event {
                ObjectEvent::ConstantControl(flags) => new_ctrl = flags,
                ObjectEvent::Rotate(rotate) => new_rotate = rotate,
            }
        }

        if new_ctrl.contains(ConstantControlFlags::FIRE_PRIMARY) &&
            snapshot.tick() >= last_fire + 20
        {
            new_fire = snapshot.tick();
            for angle in -8..8 {
                let theta = Wrapping(4096i16 * angle);
                let velocity = Affine2dH::rotate_hex(theta) * Vhd(65536, 0);

                spawn.push(UnpackedCommonObject {
                    a: self_obj.a(),
                    b: self_obj.b(),
                    theta,
                    rounded_span: 0,
                    collision_group: self_obj.collision_group(),
                    object_type: TYPE_DEBUG_PROJECTILE,
                    extended_data: 0,
                    data_dst_size: 0,
                    vax4: (velocity.a() * 4 / 100) as i16,
                    fa: 255,
                    aax16: 0,
                    vbx4: (velocity.b() * 4 / 100) as i16,
                    fb: 255,
                    abx16: 0,
                    vthetax4: 0,
                    ftheta: 255,
                    athetax16: 0,
                    wakeup_counter: 0,
                    wakeup_increment: 4,
                    id: pipeline.alloc_id(),
                }.pack());
            }
        }

        let accel = if new_ctrl.contains(ConstantControlFlags::ACCELERATE) {
            Vhd(0, 255)
        } else {
            Vhd(0, 0)
        };

        let accel = Affine2dH::rotate_hex(self_obj.theta()) * accel;
        self_obj.set_aax16(accel.a() as i8);
        self_obj.set_abx16(accel.b() as i8);
        self_obj.set_athetax16((new_rotate >> 8) as i8);

        if new_ctrl.contains(ConstantControlFlags::BRAKE) {
            self_obj.set_fa(0);
            self_obj.set_fb(0);
            self_obj.set_ftheta(0);
        } else {
            self_obj.set_fa(255);
            self_obj.set_fb(255);
            self_obj.set_ftheta(255);
        }

        if new_ctrl != old_ctrl || new_rotate != old_rotate ||
            last_fire != new_fire
        {
            let mut new_data = SmallVec::<[i32x4;32]>::from_slice(
                self_comp.as_inner());
            new_data[self_comp.as_inner().len() - 1] = i32x4::new(
                new_ctrl.bits() as i32, new_rotate as i32,
                new_fire as i32, 0);
            self_obj.set_extended_data(pipeline.alloc(&new_data));
        }

        spawn.push(self_obj);
    }
});

impl DebugPlayer {
    pub fn new(pipeline: &mut SnapshotUpdatePipeline) -> CommonObject {
        static CELLS: &'static [(i16, i16)] = &[
            (-1, -1),
            (0, -1),
            (0, 0),
            (0, 1),
            (1, -2),
        ];

        let id = pipeline.alloc_id();
        assert!(id < 256);

        let mut data: CompositeObject<SmallVec<[i32x4;32]>> = unsafe {
            CompositeObject::build(
                UnpackedCompositeHeader::default().pack(),
                || CELLS.iter().cloned())
        };
        unsafe {
            data.as_inner_mut().push(i32x4::splat(0));
        }

        UnpackedCommonObject {
            a: 0, b: 0, theta: Wrapping(0),
            vax4: 0, vbx4: 0, vthetax4: 0,
            aax16: 0, abx16: 0, athetax16: 0,
            fa: 255, fb: 255, ftheta: 255,
            rounded_span: CommonObject::round_span(data.calc_span()),
            collision_group: id as u8,
            object_type: TYPE_DEBUG_PLAYER,
            extended_data: pipeline.alloc(data.as_inner()),
            data_dst_size: data.as_inner().len() as u8,
            wakeup_counter: 255,
            wakeup_increment: 4,
            id: id,
        }.pack()
    }
}
