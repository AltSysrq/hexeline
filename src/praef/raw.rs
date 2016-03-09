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

pub use libc::c_void;
pub use libc::c_int;
pub use libc::c_char;
pub use libc::size_t;

use super::defs::*;

pub enum SimpleContext { }
pub enum System { }
pub enum MessageBus { }

pub enum Asn1NetworkIdentifierPair { }

pub type CbDrop = extern fn (this: *mut c_void);
pub type CbObjectStep = extern fn (this: *mut c_void, id: ObjectId,
                                   context: *const SimpleContext);
pub type CbObjectRewind = extern fn (this: *mut c_void, id: ObjectId,
                                     instant: Instant);
pub type CbCreateNodeObject = extern fn (
    dst: *mut c_void,
    drop: *mut CbDrop, step: *mut CbObjectStep, rewind: *mut CbObjectRewind,
    context: *const SimpleContext, id: ObjectId) -> c_int;

pub type CbEventApply = extern fn (
    object: *mut c_void, event: *const c_void,
    object_id: ObjectId, instant: Instant,
    serno: EventSerialNumber,
    context: *const SimpleContext);
pub type CbDecodeEvent = extern fn (
    dst: *mut c_void, drop: *mut CbDrop, apply: *mut CbEventApply,
    context: *const SimpleContext,
    instant: Instant, object: ObjectId, serno: EventSerialNumber,
    data: *const c_void, size: size_t) -> c_int;

pub type CbAuthIsValid = extern fn (
    context: *const SimpleContext,
    request: *const JoinRequest) -> c_int;
pub type CbAuthGen = extern fn (
    request: *mut JoinRequest,
    context: *const SimpleContext);

pub type CbPermitObjectId = extern fn (context: *const SimpleContext,
                                       id: ObjectId) -> c_int;
pub type CbAcquireId = extern fn (context: *const SimpleContext, id: ObjectId);
pub type CbDiscoverNode = extern fn (context: *const SimpleContext,
                                     netid: *const NetIdPair,
                                     id: ObjectId);
pub type CbRemoveNode = extern fn (context: *const SimpleContext,
                                   id: ObjectId);
pub type CbJoinTreeTraversed = extern fn (context: *const SimpleContext);
pub type CbHtScanProgress = extern fn (context: *const SimpleContext,
                                       numerator: u32, denominator: u32);
pub type CbAwaitingStability = extern fn (context: *const SimpleContext,
                                          node: ObjectId, systime: Instant,
                                          committed: Instant,
                                          validated: Instant);
pub type CbInformationComplete = extern fn (context: *const SimpleContext);
pub type CbClockSynced = extern fn (context: *const SimpleContext);
pub type CbGainedGrant = extern fn (context: *const SimpleContext);
pub type CbRecvUnicast = extern fn (
    context: *const SimpleContext,
    from_node: ObjectId, instant: Instant,
    data: *const c_void, size: size_t);
pub type CbLog = extern fn (context: *const SimpleContext, msg: *const c_char);

pub type CbEventOptimism = extern fn (
    context: *const SimpleContext,
    target: ObjectId,
    when: Instant,
    serno: EventSerialNumber,
    event: *const c_void) -> u32;
pub type CbEventVote = extern fn (
    context: *const SimpleContext,
    target: ObjectId,
    when: Instant,
    serno: EventSerialNumber,
    event: *const c_void) -> c_int;

#[repr(C)]
#[derive(Copy,Clone,Debug)]
pub struct VirtualNetworkLink {
    pub base_latency: u32,
    pub variable_latency: u32,
    pub firewall_grace_period: u32,
    pub reliability: u16,
    pub duplicity: u16,
}

enum VirtualNetwork { }
enum VirtualBus { }

#[link(name = "praefectus")]
extern {
    pub fn praef_simple_new(
        userdata: *mut c_void,
        bus: *mut MessageBus, selfid: *const Asn1NetworkIdentifierPair,
        std_latency: u32, profile: SystemProfile,
        ip_version: SystemIpVersion, net_locality: SystemNetworkLocality,
        mtu: u32) -> *mut SimpleContext;
    pub fn praef_simple_delete(this: *mut SimpleContext);
    pub fn praef_simple_get_system(this: *const SimpleContext) ->
        *mut System;
    pub fn praef_simple_get_userdata(this: *const SimpleContext) ->
        *mut c_void;

    pub fn praef_simple_cb_create_node_object(
        this: *mut SimpleContext,
        callback: CbCreateNodeObject,
        object_size: size_t);

    pub fn praef_simple_cb_decode_event(
        this: *mut SimpleContext,
        callback: CbDecodeEvent,
        event_size: size_t);

    pub fn praef_simple_cb_auth(
        context: *mut SimpleContext,
        is_valid: CbAuthIsValid,
        gen: CbAuthGen);

    pub fn praef_simple_cb_permit_object_id(
        context: *mut SimpleContext,
        callback: CbPermitObjectId);
    pub fn praef_simple_cb_acquire_id(
        context: *mut SimpleContext,
        callback: CbAcquireId);
    pub fn praef_simple_cb_discover_node(
        context: *mut SimpleContext,
        callback: CbDiscoverNode);
    pub fn praef_simple_cb_remove_node(
        context: *mut SimpleContext,
        callback: CbRemoveNode);
    pub fn praef_simple_cb_join_tree_traversed(
        context: *mut SimpleContext,
        callback: CbJoinTreeTraversed);
    pub fn praef_simple_cb_ht_scan_progress(
        context: *mut SimpleContext,
        callback: CbHtScanProgress);
    pub fn praef_simple_cb_awaiting_stability_t(
        context: *mut SimpleContext,
        callback: CbAwaitingStability);
    pub fn praef_simple_cb_information_complete(
        context: *mut SimpleContext,
        callback: CbInformationComplete);
    pub fn praef_simple_cb_clock_synced(
        context: *mut SimpleContext,
        callback: CbClockSynced);
    pub fn praef_simple_cb_gained_grant(
        context: *mut SimpleContext,
        callback: CbGainedGrant);
    pub fn praef_simple_cb_recv_unicast(
        context: *mut SimpleContext,
        callback: CbRecvUnicast);
    pub fn praef_simple_cb_log(
        context: *mut SimpleContext,
        callback: CbLog);

    pub fn praef_simple_cb_event_optimism(
        context: *mut SimpleContext,
        callback: CbEventOptimism);
    pub fn praef_simple_cb_event_vote_t(
        context: *mut SimpleContext,
        callback: CbEventVote);

    pub fn praef_system_bootstrap(system: *mut System);
    pub fn praef_system_connect(system: *mut System,
                                target: *mut Asn1NetworkIdentifierPair);
    pub fn praef_system_disconnect(system: *mut System);
    pub fn praef_system_advance(system: *mut System,
                                delta: u32) -> SystemStatus;
    pub fn praef_system_add_event(system: *mut System,
                                  data: *const c_void, size: size_t) -> c_int;
    // No need to expose praef_system_vote_event(), it's covered by the
    // simplified API.
    pub fn praef_system_send_unicast(system: *mut System,
                                     target: ObjectId,
                                     data: *const c_void, size: size_t) ->
        c_int;
    pub fn praef_system_get_local_id(system: *const System) -> ObjectId;
    pub fn praef_system_get_latency_to(system: *const System,
                                   node: ObjectId) -> u32;
    pub fn praef_system_oom(system: *mut System);

    pub fn praef_system_conf_clock_obsolescence_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_clock_tolerance(
        system: *mut System, val: u32);
    pub fn praef_system_conf_commit_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_max_commit_lag(
        system: *mut System, val: u32);
    pub fn praef_system_conf_max_validated_lag(
        system: *mut System, val: u32);
    pub fn praef_system_conf_commit_lag_laxness(
        system: *mut System, val: u32);
    pub fn praef_system_conf_commit_lag_compensation(
        system: *mut System,
        numerator: u32, denominator: u32);
    pub fn praef_system_conf_public_visibility_lag(
        system: *mut System, val: u32);
    pub fn praef_system_conf_stability_wait(
        system: *mut System, val: u32);
    pub fn praef_system_conf_join_tree_query_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_accept_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_max_live_nodes(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_range_max(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_range_query_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_scan_redundancy(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_scan_concurrency(
        system: *mut System, val: u8);
    pub fn praef_system_conf_ht_max_scan_tries(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_snapshot_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_num_snapshots(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_root_query_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ht_root_query_offset(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ungranted_route_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_granted_route_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_ping_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_max_pong_silence(
        system: *mut System, val: u32);
    pub fn praef_system_conf_route_kill_delay(
        system: *mut System, val: u32);
    pub fn praef_system_conf_propose_grant_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_vote_deny_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_vote_chmod_offset(
        system: *mut System, val: u32);
    pub fn praef_system_conf_grace_period(
        system: *mut System, val: u32);
    pub fn praef_system_conf_direct_ack_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_indirect_ack_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_linear_ack_interval(
        system: *mut System, val: u32);
    pub fn praef_system_conf_linear_ack_max_xmit(
        system: *mut System, val: u32);
    pub fn praef_system_conf_max_advance_per_frame(
        system: *mut System, val: u32);
    pub fn praef_system_conf_max_event_vote_offset(
        system: *mut System, val: u32);

    pub fn praef_virtual_network_new() -> *mut VirtualNetwork;
    pub fn praef_virtual_network_delete(this: *mut VirtualNetwork);
    pub fn praef_virtual_network_create_node(
        net: *mut VirtualNetwork) -> *mut VirtualBus;
    pub fn praef_virtual_bus_link(
        a: *mut VirtualBus,
        b: *mut VirtualBus) -> *mut VirtualNetworkLink;
    pub fn praef_virtual_bus_mb(
        vb: *mut VirtualBus) -> *mut MessageBus;
    pub fn praef_virtual_bus_address(vb: *const VirtualBus) ->
        *const Asn1NetworkIdentifierPair;
    pub fn praef_virtual_bus_bw_in(vb: *const VirtualBus) -> u64;
    pub fn praef_virtual_bus_bw_out(vb: *const VirtualBus) -> u64;
    pub fn praef_virtual_network_advance(
        net: *mut VirtualNetwork, delta: u32) -> c_int;
}

// TODO: udp-message-bus (needs more help to deal with
// PraefNetworkIdentifierPair_t)
