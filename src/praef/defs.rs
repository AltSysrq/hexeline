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

pub type Instant = u32;
pub type ObjectId = u32;
pub type EventSerialNumber = u32;
pub const BOOTSTRAP_NODE : ObjectId = 1;
pub const NULL_OBJECT_ID : ObjectId = 0;

#[repr(C)]
#[derive(Copy,Clone,Debug,PartialEq,Eq)]
pub enum SystemIpVersion {
    Any, Only4, Only6
}

#[repr(C)]
#[derive(Copy,Clone,Debug,PartialEq,Eq)]
pub enum SystemNetworkLocality {
    Any, Local, Global
}

#[repr(C)]
#[derive(Copy,Clone,Debug,PartialEq,Eq,PartialOrd,Ord)]
pub enum SystemStatus {
    Ok, Anonymous, PendingGrant, Partitioned, Kicked,
    Oom, Overflow, Collision,
}

#[repr(C)]
#[derive(Copy,Clone,Debug,PartialEq,Eq)]
pub enum SystemProfile {
    Strict,
    Lax,
}

#[repr(u8)]
#[derive(Copy,Clone,Debug,PartialEq,Eq)]
pub enum IpAddressVersion {
    V4 = 4,
    V6 = 6,
}

#[repr(C)]
#[derive(Copy,Clone,Debug)]
pub struct IpAddress {
    pub version: IpAddressVersion,
    pub v4: [u8; 4],
    pub v6: [u16; 8],
}

#[repr(C)]
#[derive(Copy,Clone,Debug)]
pub struct NetId {
    pub address: IpAddress,
    pub port: u16,
}

#[repr(C)]
#[derive(Copy,Clone,Debug)]
pub struct NetIdPair {
    pub global: u8,
    pub intranet: NetId,
    pub internet: NetId,
}

#[repr(C)]
pub struct JoinRequest {
    pub public_key: [u8;32],
    pub identifier: NetIdPair,
    pub auth_size: u8,
    pub auth: [u8;58],
}
