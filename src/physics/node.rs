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

use std::collections::{BTreeMap,VecDeque};

use smallvec::SmallVec;

use praef;
use super::defs::*;
use super::particle::*;
use super::pslot::*;

#[derive(Clone,Copy,Debug,PartialEq,Eq,PartialOrd,Ord)]
struct EventId {
    instant: Chronon,
    serno: praef::EventSerialNumber,
}

#[derive(Clone,Debug)]
struct EventState<E> {
    event: E,
    effective: bool,
    // Option<PslotIndex> might seem more natural, but we need to deal with the
    // possibility of a malicious node (or a bug) creating multiple particles
    // which act as addressees of the same event, since we don't have stable
    // pslot assignments. That is, either all nodes need to agree on the one
    // particle that gets the event in such cases, or we must apply such events
    // to all applicable particles.
    applied: SmallVec<[PslotIndex;1]>,
}

/// Maintains the state of a single node. It corresponds to a praefectus
/// object, but is not one itself since it is owned elsewhere.
pub struct Node<P: Particle> {
    /// The id of this node.
    id: praef::ObjectId,
    /// The current instant of this node. Rewinding the node recedes this; it
    /// is advanced by calls to next_frame().
    now: Chronon,
    /// Counter incremented every time a particle is spawned. Used to select
    /// the starting point in `slots` to probe for a free slot.
    spawn_counter: PslotIndex,
    /// The slots in this node. The vector is expanded as needed when a spawn
    /// attempt traverses more than half of the slots without finding an empty
    /// one.
    slots: Vec<Pslot<P>>,
    /// Tracks events which have been applied to particles and whether
    /// praefectus thinks they should still be applied.
    events: BTreeMap<EventId,EventState<P::Event>>,
}

impl<P: Particle> Node<P> {
    /// Creates a new node.
    ///
    /// `mother` is the eternal, primordial pseudo-particle installed into slot
    /// 0. `t` is the initial value of `now()`.
    pub fn new(id: praef::ObjectId, mother: P, t: Chronon) -> Self {
        Node {
            id: id,
            now: t,
            spawn_counter: 0,
            slots: vec![Pslot::spawn_primordial(mother, t)],
            events: BTreeMap::new(),
        }
    }

    /// Returns the id of this node.
    pub fn id(&self) -> praef::ObjectId {
        self.id
    }

    /// Returns this node's current concept of the current time.
    pub fn now(&self) -> Chronon {
        self.now
    }

    /// Sets this node's current concept of the current time. This may only be
    /// used to reduce `now()`. This causes all events between `now()` and
    /// `then` to be marked as non-effective. Essentially, this should be a
    /// direct mapping from `praef::Object::rewind()`. If `now()` is to be
    /// adjusted for other reasons, use `set_now()`, or else events may be
    /// spuriously discarded.
    pub fn rewind(&mut self, then: Chronon) {
        use std::collections::Bound::{Included,Unbounded};

        self.set_now(then);

        for (_, state) in self.events.range_mut(
            Included(&EventId { instant: then, serno: 0 }),
            Unbounded)
        {
            state.effective = false;
        }
    }

    /// Sets this node's current conecept of the current time. This may only be
    /// used to reduce `now()`. Unlike `rewind()`, events are not marked as
    /// non-effective by this call.
    pub fn set_now(&mut self, then: Chronon) {
        debug_assert!(then <= self.now);

        self.now = then;
    }

    /// Adds an event to this node's state. If the identified event is not
    /// known, `supplier` is invoked to get the event data. In either case, the
    /// event is marked as effective so it will be applied when `t` is reached.
    pub fn add_event<F>(&mut self, t: Chronon, serno: praef::EventSerialNumber,
                        supplier: F)
    where F: FnOnce() -> P::Event {
        debug_assert!(self.now <= t);

        self.events.entry(EventId { instant: t, serno: serno })
            .or_insert_with(|| EventState {
                event: supplier(),
                effective: true,
                applied: SmallVec::new(),
            }).effective = true;
    }

    /// Called at the start of a frame.
    ///
    /// `prepare_frame()` is called on all pslots with `now()`.
    ///
    /// All events in this frame are tested. If an event is applied but not
    /// effective, the pslot to which it had been applied is damaged.
    ///
    /// If any cascading damage occurs as a result of the above, `damage` is
    /// populated with additional pslots to damage.
    pub fn prepare_frame(&mut self, damage: &mut VecDeque<PslotId>) {
        use std::collections::Bound::{Included,Excluded};

        for ps in self.slots.iter_mut() {
            ps.prepare_frame(self.now);
        }

        let mut to_remove: Vec<EventId> = vec![];
        let slots = &mut self.slots;

        // The two passes here are necessary since the first pass can cause
        // pslots to become active.
        for (id, state) in self.events.range_mut(
            Included(&EventId { instant: self.now, serno: 0 }),
            Excluded(&EventId { instant: self.now + 1, serno: 0 }))
        {
            // If the event is no longer supposed to apply to anything, damage
            // the particle to which it was applied (if any) and remove it from
            // the table.
            if !state.effective {
                to_remove.push(*id);

                for ix in state.applied.iter() {
                    let ps = &mut slots[*ix as usize];
                    ps.seek(self.now);
                    ps.damage(damage);
                }
            }
        }

        for (_, state) in self.events.range_mut(
            Included(&EventId { instant: self.now, serno: 0 }),
            Excluded(&EventId { instant: self.now + 1, serno: 0 }))
        {
            // If the particle this event applied to is currently active, we'll
            // need to reapply the event later.
            if state.applied.iter().any(
                |ix| slots[*ix as usize].is_active())
            {
                // Since there could be multiple particles entangled by an id
                // conflict, make sure to damage them all.
                for ix in state.applied.iter() {
                    let ps = &mut slots[*ix as usize];
                    ps.seek(self.now);
                    ps.damage(damage);
                }

                state.applied.clear();
            }
        }

        for id in to_remove.iter() {
            self.events.remove(id);
        }
    }

    /// Damages the pslot identified by `slot` within this node.
    ///
    /// If any cascading damage occurs, the ids of further pslots to damage are
    /// added to `damage`.
    pub fn damage(&mut self, slot: PslotIndex, damage: &mut VecDeque<PslotId>) {
        let ps = &mut self.slots[slot as usize];
        ps.seek(self.now);
        ps.damage(damage);
    }

    /// Called after collision detection to finish up the current frame.
    ///
    /// Unaplied events are applied to active pslots. All active pslots then
    /// undergo their nonanalytic advance. If any particles are spawned, they
    /// are assigned to free, active pslots. Finally, `next_frame()` is called
    /// on every active pslot, and `now()` increases by 1.
    pub fn finish_frame(&mut self, env: &P::Env) {
        use std::collections::Bound::{Included,Excluded};

        for (_, state) in self.events.range_mut(
            Included(&EventId { instant: self.now, serno: 0 }),
            Excluded(&EventId { instant: self.now + 1, serno: 0 }))
        {
            debug_assert!(state.effective);

            // Need to iterate over all slots regardless of whether the event
            // appears applied, since something new that it's addressed to
            // could have come into existence.
            for (ix,ps) in self.slots.iter_mut().enumerate() {
                // Consider applying if addressed to this slot
                if ps.is_addressee_of(&state.event) &&
                    // Ignore if already applied to this slot
                    !state.applied.iter().any(|a| (ix as PslotIndex) == *a)
                {
                    ps.apply_event(&state.event, env);
                    state.applied.push(ix as PslotIndex);
                    // Need to keep going in case there's an id conflict
                }
            }
        }

        struct ToSpawn<P> {
            child: P,
            parent: PslotIndex,
        }
        let mut to_spawn: Vec<ToSpawn<P>> = vec![];

        for (ix,ps) in self.slots.iter_mut().enumerate() {
            if ps.is_active() {
                let mut spawned = ps.advance_nonanalytic(env);
                for child in spawned.into_iter() {
                    to_spawn.push(ToSpawn { child: child,
                                            parent: ix as PslotIndex });
                }
            }
        }

        for ts in to_spawn {
            self.spawn(ts.child, ts.parent);
        }

        for ps in self.slots.iter_mut() {
            if ps.is_active() {
                ps.next_frame();
            }
        }
        self.now += 1;
    }

    fn spawn(&mut self, child: P, parent: PslotIndex) {
        let mut tries_left = self.slots.len() / 2;
        let mut child_ix: Option<PslotIndex> = None;

        // Linearly probe `slots` until we either find an empty, active slot,
        // or we've scanned half of the slots.
        while tries_left > 0 && child_ix.is_none() {
            let ix = self.spawn_counter;
            self.spawn_counter += 1;
            if (self.spawn_counter as usize) == self.slots.len() {
                self.spawn_counter = 1; // 0 always taken by mother
            }

            let ps = &mut self.slots[ix as usize];
            // Only consider active slots so we don't need to deal with damage
            // propagation.
            if ps.is_active() && !ps.is_in_use() {
                child_ix = Some(ix);
            }

            tries_left -= 1;
        }

        // If we failed to find anything, create a new slot at the end and use
        // that.
        let ix = child_ix.unwrap_or_else(|| {
            self.slots.push(Pslot::empty(self.now));
            (self.slots.len() - 1) as PslotIndex
        });

        self.slots[ix as usize].spawn(child);
        // If the parent is reset before this point, this pslot must also reset
        // since it might not get spawned.
        self.slots[parent as usize].collided_by(
            PslotId { node: self.id, pslot: ix },
            self.now + 1);
    }

    pub fn slots(&self) -> &[Pslot<P>] {
        self.slots.as_slice()
    }

    pub fn slots_mut(&mut self) -> &mut [Pslot<P>] {
        self.slots.as_mut_slice()
    }
}

#[cfg(test)]
mod test {
    use std::collections::{HashSet,VecDeque};

    use smallvec::SmallVec;

    use praef;
    use super::*;
    use super::super::defs::*;
    use super::super::particle::*;
    use super::super::pslot::*;

    struct Env;
    const ENV: &'static Env = &Env;
    #[derive(Clone,Copy,Debug,Default)]
    struct Cache;

    #[derive(Clone,Debug)]
    struct TPart {
        id: u32,
        ticks: Chronon,
        events: HashSet<u32>,
        spawn: Option<u32>,
    }

    impl TPart {
        pub fn new(id: u32) -> TPart {
            TPart {
                id: id,
                ticks: -1, // So start of first frame is 0
                events: HashSet::new(),
                spawn: None,
            }
        }

        pub fn mother() -> TPart {
            TPart::new(0)
        }
    }

    #[derive(Clone,Copy,Debug)]
    struct Event {
        target: u32,
        id: u32,
        spawn: Option<u32>,
    }

    impl Event {
        pub fn blank(target: u32, id: u32) -> Event {
            Event { target: target, id: id, spawn: None }
        }

        pub fn spawn(id: u32, spawn: u32) -> Event {
            Event { target: 0, id: id, spawn: Some(spawn) }
        }

        pub fn spawn_other(target: u32, id: u32, spawn: u32) -> Event {
            Event { target: target, id: id, spawn: Some(spawn) }
        }
    }

    impl Particle for TPart {
        type Env = Env;
        type Event = Event;
        type Cache = Cache;

        fn collision_bounds<'a>(&'a self, _cache: &'a Cache)
                                -> &'a CollisionTreeNode {
            unimplemented!();
        }

        fn is_addressee_of(&self, event: &Event) -> bool {
            self.id == event.target
        }

        fn advance_analytic(&mut self, _cache: &Cache, _env: &Env,
                            delta: Chronon) {
            self.ticks += delta;
        }

        fn max_analytic_advance(&self, _cache: &Cache, _env: &Env) -> u16 {
            65535
        }

        fn collide_with(&mut self, _self_cache: &mut Cache,
                        _other: &Self, _other_cache: &Cache,
                        _env: &Env,
                        _points: &[CollisionPoint]) -> bool {
            unimplemented!();
        }

        fn apply_event(&mut self, _cache: &Cache, event: &Event, _env: &Env) {
            assert_eq!(self.id, event.target);
            assert!(self.events.insert(event.id));
            assert!(self.spawn.is_none() || event.spawn.is_none());
            self.spawn = event.spawn;
        }

        fn advance_nonanalytic(&mut self, _cache: &mut Cache)
                               -> Option<NonanalyticTransition<Self>> {
            self.spawn.take().map(|sid| {
                NonanalyticTransition {
                    exists: true,
                    warp: false,
                    spawn: {
                        let mut v = SmallVec::new();
                        v.push(TPart::new(sid));
                        v
                    },
                }
            })
        }
    }

    fn new_node(id: praef::ObjectId, t: Chronon) -> Node<TPart> {
        Node::new(id, TPart::mother(), t)
    }

    fn find_all_with_id(node: &mut Node<TPart>, id: u32) -> Vec<&TPart> {
        node.slots_mut().iter_mut()
            .filter_map(|ps| ps.curr_particle(ENV))
            .filter(|p| id == p.id)
            .collect()
    }

    fn find_with_id(node: &mut Node<TPart>, id: u32) -> &TPart {
        node.slots_mut().iter_mut()
            .filter_map(|ps| ps.curr_particle(ENV))
            .filter(|p| id == p.id)
            .next().expect("Particle with that id not found")
    }

    fn step(node: &mut Node<TPart>, count: Chronon) {
        for _ in 0..count {
            let mut damage: VecDeque<PslotId> = VecDeque::new();
            node.prepare_frame(&mut damage);
            while let Some(id) = damage.pop_front() {
                node.damage(id.pslot, &mut damage);
            }
            node.finish_frame(ENV);
        }
    }

    #[test]
    fn simple_linear_spawn() {
        let mut node = new_node(4, 0);
        node.add_event(2, 0, || Event::spawn(0, 1));
        node.add_event(3, 1, || Event::blank(1, 1));
        step(&mut node, 4);
        assert_eq!(4, node.now());

        {
            let mother = find_with_id(&mut node, 0);
            assert_eq!(1, mother.events.len());
            assert!(mother.events.contains(&0));
            assert_eq!(4, mother.ticks);
        }

        {
            let child = find_with_id(&mut node, 1);
            assert_eq!(1, child.events.len());
            assert!(child.events.contains(&1));
            assert_eq!(1, child.ticks);
        }
    }

    #[test]
    fn event_application_and_redaction() {
        let mut node = new_node(4, 0);
        node.add_event(2, 0, || Event::blank(0, 0));
        node.add_event(2, 1, || Event::blank(0, 1));
        step(&mut node, 4);

        {
            let p = find_with_id(&mut node, 0);
            assert_eq!(2, p.events.len());
            assert!(p.events.contains(&0));
            assert!(p.events.contains(&1));
            assert_eq!(4, p.ticks);
        }

        // Event 1 is redacted
        node.rewind(2);
        node.add_event(2, 0, || Event::blank(0, 0));
        step(&mut node, 2);
        assert_eq!(4, node.now());

        {
            let p = find_with_id(&mut node, 0);
            assert_eq!(1, p.events.len());
            assert!(p.events.contains(&0));
            assert_eq!(4, p.ticks);
        }
    }

    #[test]
    fn dependent_spawn_redaction() {
        let mut node = new_node(4, 0);
        // Particle 0 spawns 1 and 2 at t=1 and t=2.
        // Particle 1 spawns 3 at t=2.
        node.add_event(1, 0, || Event::spawn(0, 1));
        node.add_event(2, 1, || Event::spawn(1, 2));
        node.add_event(2, 2, || Event::spawn_other(1, 2, 3));
        step(&mut node, 4);

        assert_eq!(1, find_all_with_id(&mut node, 0).len());
        assert_eq!(1, find_all_with_id(&mut node, 1).len());
        assert_eq!(1, find_all_with_id(&mut node, 2).len());
        assert_eq!(1, find_all_with_id(&mut node, 3).len());

        // Spawning of 1 redacted. Neither 1 nor 3 should come into existence,
        // even though the event that spawned 3 is still in effect.
        node.rewind(1);
        node.add_event(1, 1, || Event::spawn(1, 2));
        node.add_event(2, 2, || Event::spawn_other(1, 2, 3));
        step(&mut node, 3);

        assert_eq!(1, find_all_with_id(&mut node, 0).len());
        assert_eq!(0, find_all_with_id(&mut node, 1).len());
        assert_eq!(1, find_all_with_id(&mut node, 2).len());
        assert_eq!(0, find_all_with_id(&mut node, 3).len());
    }
}
