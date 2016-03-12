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

#![allow(dead_code)]

use std::marker::PhantomData;

use super::raw;
use super::traits::*;

pub use super::raw::VirtualNetworkLink;

pub struct VirtualNetwork {
    vnet: *mut raw::VirtualNetwork,
}

pub struct VirtualBus<'a> {
    vbus: *mut raw::VirtualBus,
    _life: PhantomData<&'a raw::VirtualBus>,
}

impl Drop for VirtualNetwork {
    fn drop(&mut self) {
        unsafe {
            raw::praef_virtual_network_delete(self.vnet);
        }
    }
}

impl VirtualNetwork {
    pub fn new() -> Option<VirtualNetwork> {
        let vnet = unsafe {
            raw::praef_virtual_network_new()
        };

        if vnet.is_null() {
            return None;
        }

        Some(VirtualNetwork { vnet: vnet })
    }

    pub fn create_node<'a>(&'a self) -> Option<VirtualBus<'a>> {
        let vbus = unsafe {
            raw::praef_virtual_network_create_node(self.vnet)
        };

        if vbus.is_null() {
            return None;
        }

        Some(VirtualBus { vbus: vbus, _life: PhantomData })
    }

    pub fn advance(&mut self, delta: u32) -> bool {
        0 != unsafe {
            raw::praef_virtual_network_advance(self.vnet, delta)
        }
    }
}

impl<'a> VirtualBus<'a> {
    pub fn link<'b,'c>(&'c mut self, other: &'c mut VirtualBus<'b>)
                       -> Option<&'c mut VirtualNetworkLink> {
        unsafe {
            let ptr = raw::praef_virtual_bus_link(
                self.vbus, other.vbus);

            if ptr.is_null() {
                None
            } else {
                Some(&mut*ptr)
            }
        }
    }

    pub fn bw_in(&self) -> u64 {
        unsafe {
            raw::praef_virtual_bus_bw_in(self.vbus)
        }
    }

    pub fn bw_out(&self) -> u64 {
        unsafe {
            raw::praef_virtual_bus_bw_out(self.vbus)
        }
    }
}

unsafe impl<'a> MessageBus for VirtualBus<'a> {
    fn to_raw_message_bus(&self) -> *mut raw::MessageBus {
        unsafe {
            raw::praef_virtual_bus_mb(self.vbus)
        }
    }

    fn self_netid(&self) -> *const raw::Asn1NetworkIdentifierPair {
        unsafe {
            raw::praef_virtual_bus_address(self.vbus)
        }
    }
}
