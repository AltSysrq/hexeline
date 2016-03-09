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

use super::defs::*;
use super::raw;
use super::traits::*;

pub struct System<'a, S: Application + 'a> {
    context: *mut raw::SimpleContext,
    state: &'a mut S,
    connected: bool,
}

impl <'a, S: Application + 'a> System<'a, S> {
    pub fn new(application: &'a mut S,
               bus: &'a mut MessageBus,
               std_latency: u32,
               profile: SystemProfile,
               ip_version: SystemIpVersion,
               net_locality: SystemNetworkLocality,
               mtu: u32) -> Self {
        unsafe {
            let context = raw::praef_simple_new(
                application as *mut S as *mut raw::c_void,
                bus.to_raw_message_bus(), bus.get_self_netid(),
                std_latency, profile, ip_version, net_locality, mtu);
            if context.is_null() {
                panic!("Failed to allocate simplified context");
            }

            System {
                state: application,
                context: context,
                connected: false,
            }
        }
    }

    fn bootstrap(&mut self) {
        if self.connected {
            panic!("Context already connected");
        }

        self.connected = true;
        unsafe {
            raw::praef_system_bootstrap(self.sys());
        }
    }

    fn disconnect(&mut self) {
        if !self.connected {
            panic!("Context not yet connected");
        }
        unsafe {
            raw::praef_system_disconnect(self.sys());
        }
    }

    fn sys(&self) -> *mut raw::System {
        unsafe {
            raw::praef_simple_get_system(self.context)
        }
    }
}

impl<'a, S: Application + 'a> Drop for System<'a, S> {
    fn drop(&mut self) {
        unsafe {
            raw::praef_simple_delete(self.context);
        }
    }
}
