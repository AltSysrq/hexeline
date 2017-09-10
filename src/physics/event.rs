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

bitflags! {
    /// Flags for "constant controls".
    ///
    /// "Constant controls" are digital controls which take effect over time
    /// (as opposed to having an effect exactly once).
    pub struct ConstantControlFlags : u32 {
        // Movement
        const ACCELERATE        = 1 << 0;
        const BRAKE             = 1 << 1;
        const STRAFE_LEFT       = 1 << 2;
        const STRAFE_RIGHT      = 1 << 3;
        const OVERDRIVE         = 1 << 4;

        // Weapons
        const FIRE_PRIMARY      = 1 << 10;
        const FIRE_SECONDARY    = 1 << 11;
        const RELOAD            = 1 << 12;

        // Misc
        const ACCESSORIES       = 1 << 20;
    }
}

/// Events which apply to a particular object.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ObjectEvent {
    /// Reset the constant controls on the object.
    ConstantControl(ConstantControlFlags),
    /// Reset the rotational control on the object. -32768 is the maximum
    /// counter-clockwise control, 32767 the maximum clockwise control.
    Rotate(i16),
}
