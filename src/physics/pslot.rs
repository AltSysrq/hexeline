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

use std::collections::hash_map::HashMap;
use std::collections::vec_deque::VecDeque;

use praef;
use super::defs::*;
use super::particle::*;

struct TapeElement<T: Default> {
    t: Chronon,
    v: T,
}

impl<T: Default> Default for TapeElement<T> {
    fn default() -> Self {
        TapeElement {
            t: -1,
            v: Default::default()
        }
    }
}

/// Stores a sequence of items associated with a timestamp. At all times, a
/// tape has a "current" item, which is the item with the largest timestamp
/// less than the current time. Note that this means that items represent
/// values at the *end* of each instant.
struct Tape<T: Default> {
    values: Vec<TapeElement<T>>,
    ix: usize,
}

impl<T: Default> Default for Tape<T> {
    fn default() -> Self {
        Tape {
            ix: 0,
            values: vec![Default::default()],
        }
    }
}

impl<T: Default> Tape<T> {
    fn seek(&mut self, target: Chronon) {
        while self.ix > 0 && self.values[self.ix].t >= target {
            self.ix -= 1;
        }
        while self.ix + 1 < self.values.len() &&
            self.values[self.ix + 1].t < target
        {
            self.ix += 1;
        }
    }

    fn curr(&self) -> &T {
        &self.values[self.ix].v
    }

    fn curr_t(&self) -> Chronon {
        self.values[self.ix].t
    }

    fn next(&self) -> Option<&T> {
        if self.ix + 1 < self.values.len() {
            Some(&self.values[self.ix + 1].v)
        } else {
            None
        }
    }

    fn next_t(&self) -> Option<Chronon> {
        if self.ix + 1 < self.values.len() {
            Some(self.values[self.ix + 1].t)
        } else {
            None
        }
    }

    fn push(&mut self, t: Chronon, v: T) {
        debug_assert!(self.ix == self.values.len() - 1);
        self.values.push(TapeElement { t: t, v: v });
        self.ix += 1;
    }

    fn drop_future(&mut self) {
        self.values.truncate(self.ix + 1)
    }
}

#[derive(Copy,Clone,PartialEq,Eq,Hash,Debug)]
pub struct PslotId {
    pub node: praef::ObjectId,
    pub pslot: u32,
}

impl PslotId {
    fn new(node: praef::ObjectId, pslot: u32) -> PslotId {
        PslotId { node: node, pslot: pslot }
    }
}

struct PslotCurr<P: Particle> {
    /// The current state of the particle. `None` if the current state has not
    /// yet been computed.
    value: Option<P>,
    /// The cache for this particle.
    cache: P::Cache,
}

// Can't derive because that adds `P: Default` for some reason.
impl<P: Particle> Default for PslotCurr<P> {
    fn default() -> Self {
        PslotCurr {
            value: Default::default(),
            cache: Default::default(),
        }
    }
}

/// Manages the historical state of a set of particles in a single identifier
/// slot.
///
/// See the module documentation of `particle` for more details on the notional
/// functionality of pslots.
///
/// While the `Pslot` abstracts away the actual management of the data, calling
/// code must still be aware of the state machine each `Pslot` implements.
///
/// Most importantly is whether the `Pslot` is active or passive (see
/// `is_active()`). Many functions only work on active or passive `Pslot`s. A
/// passive `Pslot` is primarily manipulated with `Pslot::seek()`,
/// `Pslot::curr_particle()`, and `Pslot::approx_particle()` to read its state
/// at particular times. Active `Pslot`s are stepped through
/// `Pslot::collide_with()`, `Pslot::apply_event()`,
/// `Pslot::advance_nonanalytic()`, `Pslot::spawn()`, and
/// `Pslot::next_frame()`, in that order.
///
/// An active `Pslot` is said to be in a "transitional" state when it has
/// pending nonanalytic mutations that have yet to be saved to the state list.
/// A transitional `Pslot` cannot be `seek()`ed or otherwise manipulated in a
/// way that would lose this state; additionally, reading the current state
/// returns the transitional state rather than a pure "start-of-frame" state.
///
/// The "start-of-frame state" of a particle is the result of applying analytic
/// updates to the last non-analytic transition of that particle, up to and
/// including that frame.
pub struct Pslot<P: Particle> {
    /// The known non-analytic transitions of the pslot, past and future.
    states: Tape<Option<P>>,
    /// Warp points of the pslot, past and future.
    warps: Tape<()>,
    /// Map from observers which have collided with this pslot to the maximum
    /// time that has been seen from such a collision.
    ///
    /// There is no bookkeeping to account for rewinding time behind each
    /// value. The maximum value is always kept, even if there ends up being no
    /// collision at that point.
    max_collision_subject: HashMap<PslotId,Chronon>,
    /// The current state of the particle in this slot, or `None` if there is
    /// no particle here right now.
    curr: Option<PslotCurr<P>>,
    /// The current time of the state in this pslot.
    now: Chronon,
    /// Whether any non-analytic updates have been applied to the current
    /// particle state during the current chronon.
    nonanalytic: bool,
    /// Whether any warps have been applied to the current particle during the
    /// current chronon.
    warped: bool,
    /// The current deterministic limit, exclusive. Specifically, this refers
    /// to the first chronon at which this pslot is active. This is always
    /// greater than the greatest value of `t` in `states`.
    det_limit: Chronon,
}

// Due to https://github.com/rust-lang/rust/issues/26925 again, need to write
// Default manually.
impl<P: Particle> Default for Pslot<P> {
    fn default() -> Self {
        Pslot {
            states: Default::default(),
            warps: Default::default(),
            max_collision_subject: Default::default(),
            curr: Default::default(),
            now: 0,
            nonanalytic: false,
            warped: false,
            det_limit: 0,
        }
    }
}

impl<P: Particle> Pslot<P> {
    /// Returns whether this `Pslot` is currently active.
    pub fn is_active(&self) -> bool {
        self.now >= self.det_limit
    }

    /// Returns the current time represented by this `Pslot`.
    pub fn now(&self) -> Chronon {
        self.now
    }

    /// Returns whether this `Pslot` is currently in a transitional state.
    pub fn is_transitional(&self) -> bool {
        self.nonanalytic
    }

    /// Changes the current time of state represented by this `Pslot`.
    ///
    /// The time must be within the passive range, or else be the first
    /// upcoming active frame.
    ///
    /// This may not be invoked on active `Pslot`s in transitional states.
    pub fn seek(&mut self, t: Chronon) {
        debug_assert!(t <= self.det_limit && t >= 0);
        debug_assert!(!self.nonanalytic);

        if t != self.now {
            self.now = t;
            self.states.seek(t);
            self.warps.seek(t);
            // Ensure curr represents whether the particle exists, and mark
            // the state of the particle as non-computed.
            match *self.states.curr() {
                None => {
                    self.curr = None;
                }
                Some(_) => {
                    if self.curr.is_none() {
                        self.curr = Some(Default::default());
                    } else {
                        self.curr.as_mut().unwrap().value = None;
                    }
                }
            }
        }
    }

    /// Makes this `Pslot` active if it wasn't already.
    ///
    /// The current now of the `Pslot` becomes the new deterministc limit.
    ///
    /// `PslotId`s corresponding to other `Pslot`s that must be damaged as a
    /// result of rolling back to this time are appended to `damage_others`.
    ///
    /// This call has no effect if this `Pslot` is already active.
    pub fn damage(&mut self, damage_others: &mut VecDeque<PslotId>) {
        if !self.is_active() {
            self.det_limit = self.now;
            self.states.drop_future();
            self.warps.drop_future();

            for (id, t) in &self.max_collision_subject {
                if *t >= self.now {
                    damage_others.push_back(*id);
                }
            }
        }
    }

    fn ensure_curr_state(&mut self, env: &P::Env) {
        if let Some(c) = self.curr.as_mut() {
            if c.value.is_none() {
                let mut p = self.states.curr().clone().unwrap();
                p.advance_analytic(&c.cache, env,
                                   self.now - self.states.curr_t());
                c.value = Some(p);
            }
        }
    }

    fn map_particle<'a,T,F>(&'a mut self, env: &P::Env, f: F) -> Option<T>
    where F: FnOnce(&'a mut P, &'a mut P::Cache) -> T {
        self.ensure_curr_state(env);
        self.map_particle_ready_mut(f)
    }

    fn map_particle_ready<'a,T,F>(&'a self, f: F) -> Option<T>
    where F: FnOnce(&'a P, &'a P::Cache) -> T {
        self.curr.as_ref().map(|c| f(c.value.as_ref().unwrap(), &c.cache))
    }

    fn map_particle_ready_mut<'a,T,F>(&'a mut self, f: F) -> Option<T>
    where F: FnOnce(&'a mut P, &'a mut P::Cache) -> T {
        self.curr.as_mut().map(|c| f(c.value.as_mut().unwrap(), &mut c.cache))
    }

    /// Returns the current state of the contained particle, computing it if
    /// necessary.
    ///
    /// This may be used on both active and passive particles. For active
    /// particles in a transitional state, the resulting particle state
    /// reflects the partial transition applied so far.
    ///
    /// Returns None if there is no particle at this time.
    pub fn curr_particle(&mut self, env: &P::Env) -> Option<&P> {
        self.map_particle(env, |p, _| &*p)
    }

    /// Returns an "approximate" state of the contained particle.
    ///
    /// The returned state will reflect the exact start-of-frame state of the
    /// particle at some instant at or before `not_after`, but where no warps
    /// occur between that time and `not_after`. The actual time of the
    /// returned state is returned as the first member of the tuple.
    ///
    /// This may only be called on passive particles, and may not be used to
    /// query for states beyond their deterministic limit (since it is unknown
    /// whether warps may occur beyond that point).
    ///
    /// This call may invoke `seek()` implicitly if necessary to fullfil the
    /// request.
    pub fn approx_particle(&mut self, not_after: Chronon)
                           -> (Chronon,Option<&P>) {
        debug_assert!(!self.is_active());
        debug_assert!(not_after <= self.det_limit);

        if self.now > not_after ||
            self.warps.next_t().map_or(false, |t| t < not_after)
        {
            self.seek(not_after);
        }

        if let Some(PslotCurr { value: Some(ref val), ..}) = self.curr {
            (self.now, Some(val))
        } else {
            (self.states.curr_t(), self.states.curr().as_ref())
        }
    }

    /// Returns whether the given event is addressed to this `Pslot`, based on
    /// whether the `Pslot` is in use and whether it is addressed to the
    /// contained particle.
    pub fn is_addressee_of(&self, event: &P::Event) -> bool {
        self.states.curr().as_ref().map_or(false, |p| p.is_addressee_of(event))
    }

    /// Like `Particle::max_analytic_advance()`.
    pub fn max_analytic_advance(&mut self, env: &P::Env) -> u16 {
        self.map_particle(env, |p, cache| p.max_analytic_advance(cache, env))
            .unwrap_or(65535)
    }

    /// Returns a lower bound on when the `Pslot` may next warp.
    ///
    /// Warp times specifically indicate the chronon at which the warp occurs;
    /// thus a warp won't be visible in the start-of-frame state of the pslot
    /// until the first chronon _after_ the returned value.
    pub fn next_warp(&self) -> Chronon {
        self.warps.next_t().unwrap_or(self.det_limit)
    }

    /// Returns whether this `Pslot` is in use at its current time.
    pub fn is_in_use(&self) -> bool {
        self.curr.is_some()
    }

    /// Applies an event to the contained particle.
    ///
    /// This must only be called if `is_addressee_of()` is true, and may only
    /// be called for active particles.
    ///
    /// This puts the particle into a transitional state.
    pub fn apply_event(&mut self, event: &P::Event, env: &P::Env) {
        debug_assert!(self.is_active());
        debug_assert!(self.is_addressee_of(event));
        self.map_particle(env, |p, cache| p.apply_event(cache, event, env));
        self.nonanalytic = true;
    }

    /// Updates the particle in this `Pslot` with a collision from the particle
    /// in `other`.
    ///
    /// Both `Pslot`s must be in-use and active, and must have had their state
    /// for this frame already calculated.
    ///
    /// This may put `self` into a transitional state.
    pub fn collide_with(&mut self, other: &Self, env: &P::Env,
                        points: &[CollisionPoint]) {
        debug_assert!(self.curr.is_some());
        debug_assert!(other.curr.is_some());
        debug_assert!(self.is_active());
        debug_assert!(other.is_active());

        self.nonanalytic |= other.map_particle_ready(|other_p, other_cache| {
            self.map_particle_ready_mut(|self_p, self_cache| {
                self_p.collide_with(self_cache, other_p, other_cache,
                                    env, points)
            }).unwrap_or(false)
        }).unwrap_or(false);
    }

    /// Records that `other` collided with `self` at instant `t`, so that
    /// `damage()` can propagate to it if needed.
    pub fn collided_by(&mut self, other: PslotId, t: Chronon) {
        let val = self.max_collision_subject.entry(other).or_insert(t);
        if *val < t {
            *val = t;
        }
    }

    /// Performs the non-analytic step for this `Pslot`.
    ///
    /// If this step spawns new particles, they are added to `spawn`.
    ///
    /// This may only be called for active particles, and may put them in a
    /// transitional state.
    pub fn advance_nonanalytic(&mut self, env: &P::Env, spawn: &mut Vec<P>) {
        debug_assert!(self.is_active());

        let opt_result = self.map_particle(
            env, |p, cache| p.advance_nonanalytic(cache)).unwrap_or(None);
        if let Some(mut result) = opt_result {
            self.nonanalytic = true;
            self.warped |= result.warp;
            if !result.exists {
                self.curr = None;
                self.warped = true;
            }
            for child in result.spawn.into_iter() {
                spawn.push(child);
            }
        }
    }

    /// Spawns a particle within this pslot.
    ///
    /// The pslot must be active, and must not currently be in use. This call
    /// puts it into a transitional state.
    pub fn spawn(&mut self, particle: P) {
        debug_assert!(self.is_active());
        debug_assert!(!self.is_in_use());
        self.curr = Some(PslotCurr {
            value: Some(particle),
            cache: Default::default(),
        });
        self.warped = true;
        self.nonanalytic = true;
    }

    /// Finalises the changes to this `Pslot` for this frame, and advances to
    /// the next frame.
    ///
    /// The current instant for the `Pslot` is incremented, and its
    /// deterministic limit advances to match. If the `Pslot` was in a
    /// transitional state, the state is saved and the `Pslot` returns to a
    /// consistent state.
    ///
    /// This may only be called on active `Pslot`s.
    pub fn next_frame(&mut self) {
        debug_assert!(self.is_active());
        if self.nonanalytic {
            self.states.push(self.now, self.curr.as_mut().map(
                |c| c.value.take().unwrap()));
            self.nonanalytic = false;
        }
        if self.warped {
            self.warps.push(self.now, ());
            self.warped = false;
        }

        if let Some(curr) = self.curr.as_mut() {
            // Need to recompute from last pushed non-analytic transition.
            curr.value = None;
        }

        self.now += 1;
        self.det_limit = self.now;
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use super::super::defs::*;
    use super::super::particle::*;
    use std::collections::vec_deque::VecDeque;
    use smallvec::SmallVec;

    const INIT_CD: i32 = 4;

    struct Env;
    const ENV: &'static Env = &Env;
    #[derive(Clone,Copy,Debug,Default)]
    struct Cache;

    #[derive(Clone,Copy,Debug)]
    struct Event {
        to: i32,
        amount: i32,
        spawn_incr: i32,
    }

    impl Event {
        fn new(to: i32, amount: i32, spawn_incr: i32) -> Event {
            Event { to: to, amount: amount, spawn_incr: spawn_incr }
        }
    }

    #[derive(Clone,Copy,Debug)]
    struct TPart {
        id: i32,
        sum: i32,
        countdown: i32,
        collision_sum: i32,
        spawn: i32,
    }

    impl TPart {
        fn new(id: i32) -> TPart {
            TPart {
                id: id,
                countdown: INIT_CD,
                sum: 0,
                collision_sum: 0,
                spawn: 0,
            }
        }
    }

    impl Particle for TPart {
        type Env = Env;
        type Event = Event;
        type Cache = Cache;

        fn collision_bounds<'a>(&'a self, _cache: &'a Cache,)
                                -> &'a CollisionTreeNode {
            unimplemented!();
        }

        fn is_addressee_of(&self, event: &Event) -> bool {
            self.id == event.to
        }

        fn advance_analytic(&mut self, _cache: &Cache, _env: &Env,
                            delta: Chronon) {
            self.countdown -= delta as i32;
        }

        fn max_analytic_advance(&self, _cache: &Cache, _env: &Env) -> u16 {
            self.countdown as u16
        }

        fn collide_with(&mut self, _self_cache: &mut Cache,
                        other: &Self, _other_cache: &Cache,
                        _env: &Env,
                        _points: &[CollisionPoint]) -> bool {
            self.collision_sum += other.sum;
            other.sum != 0
        }

        fn apply_event(&mut self, _cache: &Cache, event: &Event, _env: &Env) {
            self.sum += event.amount;
            self.spawn += event.spawn_incr;
        }

        fn advance_nonanalytic(&mut self, _cache: &mut Cache)
                               -> Option<NonanalyticTransition<TPart>> {
            assert!(self.countdown >= 0);
            if 0 == self.countdown {
                self.sum += 1;
                self.countdown = INIT_CD;
                let mut ret = NonanalyticTransition {
                    exists: self.sum <= self.collision_sum,
                    warp: false,
                    spawn: SmallVec::new(),
                };
                for n in 0..self.spawn {
                    ret.spawn.push(TPart {
                        id: self.id + n + 1,
                        sum: 0,
                        spawn: 0,
                        .. *self
                    });
                };
                Some(ret)
            } else {
                None
            }
        }
    }

    #[test]
    fn simple_lifecycle() {
        let mut ps: Pslot<TPart> = Default::default();
        let mut spawn_vec: Vec<TPart> = Vec::new();

        assert_eq!(0, ps.now());
        assert!(!ps.is_transitional());
        assert!(ps.is_active());
        assert!(!ps.is_in_use());
        assert_eq!(65535, ps.max_analytic_advance(ENV));
        assert!(ps.curr_particle(ENV).is_none());

        ps.advance_nonanalytic(ENV, &mut spawn_vec);
        assert!(spawn_vec.is_empty());
        ps.spawn(TPart::new(0));
        assert!(ps.is_transitional());
        assert_eq!(INIT_CD, ps.curr_particle(ENV).unwrap().countdown);
        assert!(ps.is_in_use());
        ps.next_frame();

        for n in 1..INIT_CD+1 {
            assert_eq!(n, ps.now());
            assert!(!ps.is_transitional());
            assert!(ps.is_active());
            assert!(ps.is_in_use());
            assert_eq!(INIT_CD - n, ps.curr_particle(ENV).unwrap().countdown);

            ps.advance_nonanalytic(ENV, &mut spawn_vec);
            assert!(spawn_vec.is_empty());
            ps.next_frame();
        }

        assert_eq!(INIT_CD+1, ps.now());
        assert!(!ps.is_transitional());
        assert!(ps.is_active());
        assert!(!ps.is_in_use());
        assert!(ps.curr_particle(ENV).is_none());
        ps.advance_nonanalytic(ENV, &mut spawn_vec);
        assert!(spawn_vec.is_empty());
        ps.next_frame();

        assert_eq!(INIT_CD+2, ps.now());
        assert!(!ps.is_transitional());
        assert!(ps.is_active());
        assert!(!ps.is_in_use());
        assert!(ps.curr_particle(ENV).is_none());
    }

    fn step_cycles(ps: &mut Pslot<TPart>, n: u32) {
        let mut spawn_vec: Vec<TPart> = Vec::new();

        for _ in 0..n {
            ps.advance_nonanalytic(ENV, &mut spawn_vec);
            ps.next_frame();
        }

        assert!(spawn_vec.is_empty());
    }

    #[test]
    fn analytic_seek() {
        let mut ps: Pslot<TPart> = Default::default();

        step_cycles(&mut ps, 5);
        assert_eq!(5, ps.now());
        assert!(ps.is_active());
        assert!(!ps.is_in_use());
        ps.seek(3);
        assert_eq!(3, ps.now());
        assert!(!ps.is_active());
        assert!(!ps.is_in_use());
        ps.seek(2);
        assert_eq!(2, ps.now());
        assert!(!ps.is_active());
        assert!(!ps.is_in_use());
        ps.seek(4);
        assert_eq!(4, ps.now());
        assert!(!ps.is_active());
        assert!(!ps.is_in_use());
        ps.seek(5);
        assert_eq!(5, ps.now());
        assert!(ps.is_active());
        assert!(!ps.is_in_use());

        ps.spawn(TPart::new(0));
        step_cycles(&mut ps, 3);
        assert_eq!(8, ps.now());
        assert!(ps.is_in_use());
        assert_eq!(INIT_CD - 3, ps.curr_particle(ENV).unwrap().countdown);
        ps.seek(7);
        assert_eq!(7, ps.now());
        assert!(ps.is_in_use());
        assert!(!ps.is_active());
        assert_eq!(INIT_CD - 2, ps.curr_particle(ENV).unwrap().countdown);
        ps.seek(6);
        assert_eq!(6, ps.now());
        assert!(ps.is_in_use());
        assert!(!ps.is_active());
        assert_eq!(INIT_CD - 1, ps.curr_particle(ENV).unwrap().countdown);
        ps.seek(7);
        assert_eq!(7, ps.now());
        assert!(ps.is_in_use());
        assert!(!ps.is_active());
        assert_eq!(INIT_CD - 2, ps.curr_particle(ENV).unwrap().countdown);
        ps.seek(8);
        assert_eq!(8, ps.now());
        assert!(ps.is_in_use());
        assert!(ps.is_active());
        assert_eq!(INIT_CD - 3, ps.curr_particle(ENV).unwrap().countdown);

        step_cycles(&mut ps, 7);
        assert_eq!(15, ps.now());
        assert!(!ps.is_in_use());
        ps.seek(13);
        assert_eq!(13, ps.now());
        assert!(!ps.is_in_use());
        assert!(!ps.is_active());
        ps.seek(15);
        assert_eq!(15, ps.now());
        assert!(ps.is_active());
        assert!(!ps.is_in_use());
    }

    #[test]
    fn nonanalytic_seek() {
        let mut ps: Pslot<TPart> = Default::default();

        step_cycles(&mut ps, 4);
        ps.spawn(TPart::new(0));
        ps.next_frame();
        step_cycles(&mut ps, 2);
        ps.apply_event(&Event::new(0, -2, 0), ENV);
        ps.next_frame();
        assert_eq!(8, ps.now());
        step_cycles(&mut ps, 12);
        assert_eq!(20, ps.now());
        assert!(!ps.is_in_use());

        ps.seek(2);
        assert!(!ps.is_active());
        assert!(!ps.is_in_use());
        assert_eq!(4, ps.next_warp());

        ps.seek(5);
        assert!(!ps.is_active());
        assert!(ps.is_in_use());
        assert_eq!(INIT_CD - 1, ps.curr_particle(ENV).unwrap().countdown);
        assert_eq!(0, ps.curr_particle(ENV).unwrap().sum);

        ps.seek(6);
        assert!(!ps.is_active());
        assert!(ps.is_in_use());
        assert_eq!(INIT_CD - 2, ps.curr_particle(ENV).unwrap().countdown);
        assert_eq!(0, ps.curr_particle(ENV).unwrap().sum);

        ps.seek(9);
        assert!(!ps.is_active());
        assert!(ps.is_in_use());
        assert_eq!(INIT_CD - 1, ps.curr_particle(ENV).unwrap().countdown);
        assert_eq!(-1, ps.curr_particle(ENV).unwrap().sum);
        assert_eq!(16, ps.next_warp());

        ps.seek(16);
        assert!(!ps.is_active());
        assert!(ps.is_in_use());
        assert_eq!(0, ps.curr_particle(ENV).unwrap().countdown);
        assert_eq!(0, ps.curr_particle(ENV).unwrap().sum);
        assert_eq!(16, ps.next_warp());

        ps.seek(17);
        assert!(!ps.is_active());
        assert!(!ps.is_in_use());
        assert_eq!(20, ps.next_warp());

        ps.seek(20);
        assert!(ps.is_active());
        assert!(!ps.is_in_use());
        assert_eq!(20, ps.next_warp());
    }

    #[test]
    fn damage() {
        let mut ps: Pslot<TPart> = Default::default();
        let id1 = PslotId::new(0, 0);
        let id2 = PslotId::new(1, 1);
        let mut damaged: VecDeque<PslotId> = VecDeque::new();

        step_cycles(&mut ps, 10);
        ps.collided_by(id1, 4);
        ps.collided_by(id2, 3);
        ps.collided_by(id2, 5);
        ps.collided_by(id2, 2);

        ps.seek(5);
        assert!(!ps.is_active());
        ps.damage(&mut damaged);
        assert_eq!(1, damaged.len());
        assert_eq!(id2, damaged[0]);

        damaged.clear();
        ps.damage(&mut damaged);
        assert!(damaged.is_empty());
    }
}
