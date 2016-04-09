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
use cg;

use super::super::praef;

/// Angles are defined as a simple 32-bit signed value, where 0 corresponds to
/// zero degrees, 0x40000000 to 90 degrees counterlockwise, and so
/// forth. The use of wrapping arithmetic means that angles are always
/// normalised to [-pi,pi) radians.
pub type Angle = Wrapping<i32>;
/// Type for rotational velocity per Chronon.
pub type RotationVelocity = i32;
/// Type used for describing spatial coordinates, distances between them,
/// velocity, etc.
///
/// Coordinates are measured in microcells (uc, see CELL), which are not
/// intended to directly correspond to any real-world unit.
pub type Coord = i32;
pub type Position = cg::Vector2<Coord>;
pub type Dimension = cg::Vector2<Coord>;
pub type Speed = Coord;
pub type Velocity = cg::Vector2<Coord>;
pub type Acceleration = Coord;
/// Unit of time, equal to 1/SECOND of a real second. This is also the unit
/// used to measure Praefectus instants.
pub type Chronon = i32;

/// The "standard" small unit of measurement. This is the size of one ship cell.
pub const CELL : Coord = 1_000_000;
/// The "standard" large unit of measurement. This is about the width of the
/// screen at "average" zoom.
pub const SCREEN : Coord = CELL * 100;
/// A second expressed in chronons.
pub const SECOND : Chronon = 60;
/// The maximum speed at which anything is allowed to move. This value is chosen
/// both to reduce rewind/replay cycles on objects as past events enter the
/// system, and to ensure that chain-sized objects cannot pass through each
/// other. (Though two objects moving at the speed of light in other directions
/// potentially could here.)
pub const SPEED_OF_LIGHT : Speed = SCREEN;

pub const DEG_90 : Angle = Wrapping(16384 << 16);
pub const DEG_180 : Angle = Wrapping(-32768 << 16);
pub const DEG_270 : Angle = Wrapping(-16384 << 16);

/// Pseudo-node containing the state for the environment, ie, the game state
/// independent of any particular node. This is therefore a reserved node id.
pub const ENVIRONMENT : praef::ObjectId = 2;
