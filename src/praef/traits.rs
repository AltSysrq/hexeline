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

pub trait Object<S>: Sized {
    fn create(id: ObjectId, state: &mut S) -> Self;
    fn step(&mut self, id: ObjectId, state: &mut S);
    fn rewind(&mut self, id: ObjectId, instant: Instant);
}

#[allow(unused_variables)]
pub trait Event<S, Target: Object<S>>: Sized {
    fn decode(instant: Instant, object: ObjectId,
              serno: EventSerialNumber,
              state: &mut S,
              data: &[u8]) -> Option<Self>;

    fn apply(&self, target: &mut Target, state: &mut S,
             object_id: ObjectId, instant: Instant,
             serno: EventSerialNumber);

    fn optimism(&self, object: ObjectId, when: Instant,
                serno: EventSerialNumber,
                state: &mut S) -> u32 { 0 }
    fn should_vote(&self, object: ObjectId, when: Instant,
                   serno: EventSerialNumber,
                   state: &mut S) -> bool { true }
}

#[allow(unused_variables)]
pub trait Application: Sized {
    type O: Object<Self>;
    type E: Event<Self, Self::O>;

    fn auth_is_valid(&mut self, request: &JoinRequest) -> bool { true }
    fn auth_gen(&mut self, request: &mut JoinRequest) { }

    fn permit_object_id(&mut self, id: ObjectId) -> bool { true }
    fn acquire_id(&mut self, id: ObjectId) { }
    fn discover_node(&mut self, netid: &NetIdPair, id: ObjectId) { }
    fn remove_node(&mut self, id: ObjectId) { }
    fn join_tree_traversed(&mut self) { }
    fn ht_scan_progress(&mut self, num: u32, denom: u32) { }
    fn awaiting_stability(&mut self, node: ObjectId, systime: Instant,
                          committed: Instant, validated: Instant) { }
    fn information_complete(&mut self) { }
    fn clock_synced(&mut self) { }
    fn gained_grant(&mut self) { }
    fn recv_unicast(&mut self, from_node: ObjectId, instant: Instant,
                    data: &[u8]) { }
    fn log(&mut self, msg: &str) { }
}

pub unsafe trait MessageBus {
    fn to_raw_message_bus(&self) -> *mut raw::MessageBus;
    fn get_self_netid(&self) -> *const raw::Asn1NetworkIdentifierPair;
}
