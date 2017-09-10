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

use std::mem;

use super::common_object::CommonObject;
use super::snapshot::*;

/// Extended behaviour attached to a `CommonObject`.
///
/// In practise, all implementations are a newtype around `CommonObject`.
///
/// ## Unsafety
///
/// It is the implementor's responsibility to ensure that its arena pointer (if
/// it has one) is valid.
///
/// All implementations must be a newtype around `CommonObject`.
pub unsafe trait ExtendedObject {
    /// Perform a single-tick non-standard update on this object.
    ///
    /// This is called after `CommonObject::tick()` if `wakeup_counter` becomes
    /// zero.
    ///
    /// `collisions` reflects any collisions that happened last frame (which
    /// all have the `subject` equal to the index of this object).
    ///
    /// Collisions implicitly wake objects up. Events do _not_; an object
    /// expecting to receive events next frame must ensure that its wake-up
    /// counter is 255 at the end of this frame.
    ///
    /// `pipeline` can be used to allocate fresh memory in the arena of the
    /// snapshot being and to generate ids for newly spawned objects. It is OK
    /// to leak ids generated this way; they will be reclaimed on the next
    /// garbage collection cycle.
    ///
    /// It is the callee's responsibility to add itself to `spawn` to ensure
    /// its continued existence. If callee does not add itself to `spawn`, it
    /// will cease to exist.
    ///
    /// `snapshot` is the snapshot of the previous frame and may be examined
    /// arbitrarily.
    fn update(&self,
              pipeline: &mut SnapshotUpdatePipeline,
              spawn: &mut SpawnList,
              collisions: &[ImpendingCollision],
              snapshot: &Snapshot);
}

macro_rules! eo_deftype {
    ($name:ident $impl:tt) => {
        #[derive(Debug, Clone, Copy)]
        pub struct $name(CommonObject);
        unsafe impl ExtendedObject for $name $impl
    }
}

pub mod debug;

pub const TYPE_DEBUG_PROJECTILE: u8 = 0;
pub const TYPE_DEBUG_PLAYER: u8 = 128;

impl ExtendedObject {
    /// Cast the given `CommonObject` into the correct `ExtendedObject` type.
    pub fn from(common: &CommonObject) -> &ExtendedObject {
        macro_rules! switch {
            ($($typ:ident => $res:path,)*) => {
                match common.object_type() {
                    $($typ => {
                        let cast: &$res = unsafe {
                            mem::transmute(common)
                        };
                        cast
                    },)*
                    v => panic!("Unhandled object type {}", v),
                }
            }
        }
        switch! {
            TYPE_DEBUG_PROJECTILE => debug::DebugProjectile,
            TYPE_DEBUG_PLAYER => debug::DebugPlayer,
        }
    }

    /// Cast the given `CommonObject` into the correct `ExtendedObject` type.
    pub fn from_mut(common: &mut CommonObject) -> &mut ExtendedObject {
        macro_rules! switch {
            ($($typ:ident => $res:path,)*) => {
                match common.object_type() {
                    $($typ => {
                        let cast: &mut $res = unsafe {
                            mem::transmute(common)
                        };
                        cast
                    },)*
                    v => panic!("Unhandled object type {}", v),
                }
            }
        }
        switch! {
            TYPE_DEBUG_PROJECTILE => debug::DebugProjectile,
            TYPE_DEBUG_PLAYER => debug::DebugPlayer,
        }
    }
}
