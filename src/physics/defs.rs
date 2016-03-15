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

//! Common definitions for the physics system.

#![allow(dead_code)]

use std::num::Wrapping;

use super::super::praef;

/// Angles are defined as a simple 64-bit signed value, where 0 corresponds to
/// zero degrees, 0x4000000000000000 to 90 degrees counterlockwise, and so
/// forth. The use of wrapping arithmetic means that angles are always
/// normalised to [-pi,pi) radians.
pub type Angle = Wrapping<i64>;
/// Type for rotational velocity per Chronon.
pub type RotationVelocity = i64;
/// Type used for describing spatial coordinates, distances between them,
/// velocity, etc.
///
/// Coordinates are measured in microchains (uc, see CHAIN), which are not
/// intended to directly correspond to any real-world unit.
pub type Coord = i64;
pub type Velocity = Coord;
pub type Acceleration = Coord;
/// Unit of time, the real millisecond. This is also the unit used to measure
/// Praefectus instants. Note that this is a signed 64-bit integer so it can
/// easily be multiplied by Velocity and RotationVelocity.
pub type Chronon = i64;

/// The "standard" small unit of measurement. This is the size of one ship cell.
/// A real-world chain is 66 feet.
pub const CHAIN : Coord = 1_000_000;
/// The "standard" large unit of measurement. This is about the width of the
/// screen at "average" zoom. The definition of chain actually makes this closer
/// to a real-world nautical mile than a real-world league.
pub const LEAGUE : Coord = CHAIN * 100;
/// A second expressed in chronons.
pub const SECOND : Chronon = 1000;
/// Convenience constant for clarity.
pub const MILLISECOND : Chronon = 1;
/// The maximum speed at which anything is allowed to move. This value is chosen
/// both to reduce rewind/replay cycles on objects as past events enter the
/// system, and to ensure that chain-sized objects cannot pass through each
/// other. (Though two objects moving at the speed of light in other directions
/// potentially could here.)
///
/// This comes out to 10 leagues/sec or 36,000 knots, rather slower than the
/// speed of light in real life.
pub const SPEED_OF_LIGHT : Velocity = CHAIN / MILLISECOND;

pub const DEG_90 : Angle = Wrapping(16384 << 48);
pub const DEG_180 : Angle = Wrapping(-32768 << 48);
pub const DEG_270 : Angle = Wrapping(-16384 << 48);

/// Pseudo-node containing the state for the environment, ie, the game state
/// independent of any particular node. This is therefore a reserved node id.
pub const ENVIRONMENT : praef::ObjectId = 2;
