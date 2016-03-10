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

use std::{ptr,mem,slice,str,error,fmt};
use libc;

use super::defs::*;
use super::raw;
use super::traits::*;
use super::raw::{c_void,c_int,c_char,size_t};

pub struct System<'a, S: Application + 'a> {
    context: *mut raw::SimpleContext,
    state: &'a mut S,
    connected: bool,
}

#[derive(Copy,Clone,Debug)]
pub struct Error {
    d: &'static str,
}

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(formatter, "{}", self.d)
    }
}

impl error::Error for Error {
    fn description(&self) -> &str {
        self.d
    }
}

impl <'a, S: Application + 'a> System<'a, S> {
    pub fn new(application: &'a mut S,
               bus: &'a mut MessageBus,
               std_latency: u32,
               profile: SystemProfile,
               ip_version: SystemIpVersion,
               net_locality: SystemNetworkLocality,
               mtu: u32) -> Result<Self,Error> {
        unsafe {
            let context = raw::praef_simple_new(
                application as *mut S as *mut raw::c_void,
                bus.to_raw_message_bus(), bus.get_self_netid(),
                std_latency, profile, ip_version, net_locality, mtu);
            if context.is_null() {
                return Err(Error{ d: "Failed to allocate simplified context" });
            }

            raw::praef_simple_cb_create_node_object(
                context, Self::cb_create_node_object,
                mem::size_of::<S::O>());
            raw::praef_simple_cb_decode_event(
                context, Self::cb_decode_event,
                mem::size_of::<S::E>());
            raw::praef_simple_cb_auth(
                context, Self::cb_auth_is_valid, Self::cb_auth_gen);
            raw::praef_simple_cb_permit_object_id(
                context, Self::cb_permit_object_id);
            raw::praef_simple_cb_acquire_id(
                context, Self::cb_acquire_id);
            raw::praef_simple_cb_discover_node(
                context, Self::cb_discover_node);
            raw::praef_simple_cb_remove_node(
                context, Self::cb_remove_node);
            raw::praef_simple_cb_join_tree_traversed(
                context, Self::cb_join_tree_traversed);
            raw::praef_simple_cb_ht_scan_progress(
                context, Self::cb_ht_scan_progress);
            raw::praef_simple_cb_awaiting_stability(
                context, Self::cb_awaiting_stability);
            raw::praef_simple_cb_information_complete(
                context, Self::cb_information_complete);
            raw::praef_simple_cb_clock_synced(
                context, Self::cb_clock_synced);
            raw::praef_simple_cb_gained_grant(
                context, Self::cb_gained_grant);
            raw::praef_simple_cb_recv_unicast(
                context, Self::cb_recv_unicast);
            raw::praef_simple_cb_log(
                context, Self::cb_log);
            raw::praef_simple_cb_event_optimism(
                context, Self::cb_event_optimism);
            raw::praef_simple_cb_event_vote(
                context, Self::cb_event_vote);

            Ok(System {
                state: application,
                context: context,
                connected: false,
            })
        }
    }

    fn require_connected(&self) {
        if !self.connected {
            panic!("Context not yet connected");
        }
    }

    fn require_not_connected(&self) {
        if self.connected {
            panic!("Context already connected");
        }
    }

    pub fn bootstrap(&mut self) {
        self.require_not_connected();

        self.connected = true;
        unsafe {
            raw::praef_system_bootstrap(self.sys());
        }
    }

    pub fn disconnect(&mut self) {
        self.require_connected();
        unsafe {
            raw::praef_system_disconnect(self.sys());
        }
    }

    pub fn advance(&mut self, delta: u32) -> SystemStatus {
        self.require_connected();
        unsafe {
            raw::praef_system_advance(self.sys(), delta)
        }
    }

    pub fn add_event(&mut self, data: &[u8]) -> Result<(),Error> {
        self.require_connected();

        let success = 0 == unsafe {
            raw::praef_system_add_event(
                self.sys(), data.as_ptr() as *const c_void,
                data.len() as size_t)
        };

        if success {
            Err(Error { d: "praef_system_add_event returned failure" })
        } else {
            Ok(())
        }
    }

    pub fn send_unicast(&mut self, target: ObjectId,
                        data: &[u8]) -> Result<(),Error> {
        self.require_connected();

        let success = 0 == unsafe {
            raw::praef_system_send_unicast(
                self.sys(), target,
                data.as_ptr() as *const c_void,
                data.len() as size_t)
        };

        if success {
            Err(Error { d: "praef_system_send_unicast returned failure" })
        } else {
            Ok(())
        }
    }

    pub fn local_id(&self) -> ObjectId {
        unsafe {
            raw::praef_system_get_local_id(self.sys())
        }
    }

    pub fn latency_to(&self, to: ObjectId) -> u32 {
        unsafe {
            raw::praef_system_get_latency_to(self.sys(), to)
        }
    }

    fn sys(&self) -> *mut raw::System {
        unsafe {
            raw::praef_simple_get_system(self.context)
        }
    }

    unsafe extern fn cb_object_drop(vobject: *mut c_void) {
        Box::from_raw(vobject as *mut S::O);
    }

    unsafe extern fn cb_object_step(vobject: *mut c_void, id: ObjectId,
                                    context: *const raw::SimpleContext) {
        let object = &mut*(vobject as *mut S::O);
        object.step(id, Self::application(&context));
    }

    unsafe extern fn cb_object_rewind(vobject: *mut c_void, id: ObjectId,
                                      instant: Instant) {
        let object = &mut*(vobject as *mut S::O);
        object.rewind(id, instant);
    }

    unsafe extern fn cb_create_node_object(
        vdst: *mut c_void,
        drop: *mut raw::CbDrop, step: *mut raw::CbObjectStep,
        rewind: *mut raw::CbObjectRewind,
        context: *const raw::SimpleContext, id: ObjectId
    ) -> c_int {
        let dst = vdst as *mut S::O;
        ptr::write(dst, S::O::create(id, Self::application(&context)));
        ptr::write(drop, Self::cb_object_drop);
        ptr::write(step, Self::cb_object_step);
        ptr::write(rewind, Self::cb_object_rewind);
        1
    }

    unsafe extern fn cb_event_drop(vevent: *mut c_void) {
        Box::from_raw(vevent as *mut S::E);
    }

    unsafe extern fn cb_event_apply(
        vobject: *mut c_void, vevent: *const c_void,
        object_id: ObjectId, instant: Instant,
        serno: EventSerialNumber,
        context: *const raw::SimpleContext
    ) {
        let object = &mut*(vobject as *mut S::O);
        let event = &*(vevent as *const S::E);
        event.apply(object, Self::application(&context),
                    object_id, instant, serno);
    }

    unsafe extern fn cb_decode_event(
        vdst: *mut c_void, drop: *mut raw::CbDrop,
        apply: *mut raw::CbEventApply,
        context: *const raw::SimpleContext,
        instant: Instant, object: ObjectId, serno: EventSerialNumber,
        data: *const c_void, size: size_t
    ) -> c_int {
        let dst = vdst as *mut S::E;
        let decoded = S::E::decode(
            instant, object, serno, Self::application(&context),
            slice::from_raw_parts(data as *const u8, size as usize));

        match decoded {
            None => 0,
            Some(event) => {
                ptr::write(dst, event);
                ptr::write(drop, Self::cb_event_drop);
                ptr::write(apply, Self::cb_event_apply);
                1
            }
        }
    }

    unsafe extern fn cb_auth_is_valid(
        context: *const raw::SimpleContext,
        request: *const JoinRequest
    ) -> c_int {
        Self::application(&context).auth_is_valid(&*request) as c_int
    }

    unsafe extern fn cb_auth_gen(
        request: *mut JoinRequest,
        context: *const raw::SimpleContext
    ) {
        Self::application(&context).auth_gen(&mut*request);
    }

    unsafe extern fn cb_permit_object_id(
        context: *const raw::SimpleContext,
        id: ObjectId
    ) -> c_int {
        Self::application(&context).permit_object_id(id) as c_int
    }

    unsafe extern fn cb_acquire_id(
        context: *const raw::SimpleContext,
        id: ObjectId
    ) {
        Self::application(&context).acquire_id(id);
    }

    unsafe extern fn cb_discover_node(
        context: *const raw::SimpleContext,
        netid: *const NetIdPair,
        id: ObjectId
    ) {
        Self::application(&context).discover_node(&*netid, id);
    }

    unsafe extern fn cb_remove_node(
        context: *const raw::SimpleContext,
        id: ObjectId
    ) {
        Self::application(&context).remove_node(id);
    }

    unsafe extern fn cb_join_tree_traversed(
        context: *const raw::SimpleContext
    ) {
        Self::application(&context).join_tree_traversed();
    }

    unsafe extern fn cb_ht_scan_progress(
        context: *const raw::SimpleContext,
        numerator: u32, denominator: u32
    ) {
        Self::application(&context).ht_scan_progress(numerator, denominator);
    }

    unsafe extern fn cb_awaiting_stability(
        context: *const raw::SimpleContext,
        node: ObjectId, systime: Instant,
        committed: Instant, validated: Instant
    ) {
        Self::application(&context).awaiting_stability(
            node, systime, committed, validated);
    }

    unsafe extern fn cb_information_complete(
        context: *const raw::SimpleContext
    ) {
        Self::application(&context).information_complete();
    }

    unsafe extern fn cb_clock_synced(
        context: *const raw::SimpleContext
    ) {
        Self::application(&context).clock_synced();
    }

    unsafe extern fn cb_gained_grant(
        context: *const raw::SimpleContext
    ) {
        Self::application(&context).gained_grant();
    }

    unsafe extern fn cb_recv_unicast(
        context: *const raw::SimpleContext,
        from_node: ObjectId, instant: Instant,
        data: *const c_void, size: size_t
    ) {
        Self::application(&context).recv_unicast(
            from_node, instant,
            slice::from_raw_parts(data as *const u8, size as usize));
    }

    unsafe extern fn cb_log(
        context: *const raw::SimpleContext,
        cmsg: *const c_char
    ) {
        // Always ASCII
        let msg = str::from_utf8_unchecked(
            slice::from_raw_parts(cmsg as *const u8,
                                  libc::strlen(cmsg)));
        Self::application(&context).log(msg);
    }

    unsafe extern fn cb_event_optimism(
        context: *const raw::SimpleContext,
        target: ObjectId,
        when: Instant,
        serno: EventSerialNumber,
        vevent: *const c_void
    ) -> u32 {
        (&*(vevent as *const S::E)).optimism(
            target, when, serno, Self::application(&context))
    }

    unsafe extern fn cb_event_vote(
        context: *const raw::SimpleContext,
        target: ObjectId,
        when: Instant,
        serno: EventSerialNumber,
        vevent: *const c_void
    ) -> c_int {
        (&*(vevent as *const S::E)).should_vote(
            target, when, serno, Self::application(&context)) as c_int
    }

    unsafe fn application<'b>(context: &'b *const raw::SimpleContext) -> &'b mut S {
        let userdata = raw::praef_simple_get_userdata(*context);
        let ptr = userdata as *mut S;
        &mut*ptr
    }
}

impl<'a, S: Application + 'a> Drop for System<'a, S> {
    fn drop(&mut self) {
        unsafe {
            raw::praef_simple_delete(self.context);
        }
    }
}
