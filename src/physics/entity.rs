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

use std::collections::btree_map::BTreeMap;
use std::cmp::{Eq,Ord};

use super::*;
use praef;

/// A thunk representing the children of a node in a collision tree.
///
/// This is used for larger objects to avoid transforming their entire
/// collision tree every time fine collision must be tested.
pub trait CollisionThunk {
    /// Forces this thunk in the collision tree, producing a slice of the child
    /// nodes of the node containing this thunk.
    fn force<'a>(&'a self) -> &'a [CollisionTree<'a>];
}

/// Describes the collision boundaries of an entity.
///
/// There exist multiple types of "collisions", represented by different bits
/// in a 64-bit mask. Every tree has a "subject mask" and a "target mask";
/// given two different trees A and B which overlap spatially, A is considered
/// to collide with B if the AND of A's target mask with B's subject mask is
/// non-zero. Note that this means collisions may be non-mutual. This is
/// because "collisions" are more general, and include things such as creatures
/// looking for food (which entails the creature's sight colliding with the
/// food, but not vice-versa, since the food will not change state in
/// response).
///
/// All collision tree nodes are represented by axis-aligned bounding boxes.
///
/// Nodes in the tree may be leaves or branches. If a branch node collides with
/// something, it must be expanded and each child tested instead. If a leaf
/// node collides with something, that collision gets passed to the owning
/// entity.
pub struct CollisionTree<'a> {
    /// The bitmask of collision types of which this tree may be on the subject
    /// (or passive) end of.
    pub subject_mask: i64,
    /// The bitmask of collision types which this tree targets (ie, the active
    /// side of a collision).
    pub target_mask: i64,
    /// The centre of this bounding box. This is also the reported collision
    /// site if this is a leaf node.
    pub centre: Position,
    /// The size of the bounding box, measured from centre to each edge.
    pub size: Dimension,
    /// If Some, a thunk that can be evaluated to obtain this branch node's
    /// children. If None, this node is a leaf node.
    pub children: Option<&'a CollisionThunk>,
}

/// Describes a pairwise location of a collision between two entities.
pub struct CollisionLocation {
    /// The location of a node in the collision tree where `self` experienced
    /// the collision.
    pub self_pos: Position,
    /// The location of a node in the collision tree where `that` experienced
    /// the collision.
    pub that_pos: Position,
    /// The result of `self`'s `target_mask` ANDed with `that`'s
    /// `subject_mask`, allowing the collider to know what type(s) of
    /// collisions have occurred.
    pub mask: i64,
}

/// Describes how an entity's internal or external state changed
/// non-analytically in the process of handling an event.
pub enum EntityMutation<T> {
    /// The entity's internal state has not changed. The caller must continue
    /// using the same analytic base it has been for further updates to the
    /// object.
    Unchanged,
    /// The entity's internal state has changed. The caller must use the new
    /// state as the base for further updates.
    ///
    /// If a vector is returned, it contains a sequence of new entities which
    /// have been spawned as a result of the mutation.
    Mutated(Option<Vec<T>>),
    /// The entity has ceased to exist. This is a non-analytic transition that
    /// occurs outside the entity itself. If a vector is returned, it contains
    /// a sequence of new entities which have been spawned as a result of the
    /// mutation.
    Destroyed(Option<Vec<T>>),
}

/// Trait for the type managed by an Entity.
///
/// This essentially models a single, rigid, physical object. It defines how
/// the state of that object evolves over time, how it responds to collisions
/// and events, and fundamentally its position and basic size.
///
/// There is an important restriction on the combined behaviour of position and
/// size: No edge of the implied bounding box `position ± size` may _advance_
/// at a speed greater than `SPEED_OF_LIGHT`. There is no requirement that
/// position not move above the `SPEED_OF_LIGHT` (as long as the size recesses
/// enough such that the edge does not exceed that speed), nor is there any
/// restriction on how quickly the edges of the bounding box may _recess_.
pub trait EntityType: Clone + Sized + Eq {
    /// The type of the passive, immutable environment which the entity can
    /// observe at any time.
    type Env;
    /// The type of event that can be applied to the entity.
    type Event;

    /// Performs an "analytic" step or series of steps forward of this entity.
    ///
    /// Analytic updates are those which can be computed in constant time from
    /// a given base, regardless of the time delta and without serious loss of
    /// accuracy. For example, updating position based on velocity.
    ///
    /// Series of analytic updates are never compounded; rather, it is always
    /// computed on top of the most recent non-analytic "base" value, with an
    /// increasing time delta. Because of this, there is no requirement that
    /// ```
    ///   a.step_analytic(env, 1); a.step_analyitc(env, 1);
    /// ```
    /// produce the same result as
    /// ```
    ///   a.step_analytic(env, 2);
    /// ```
    ///
    /// Analytic steps have no way to interact with other entities; the only
    /// external information they can see is the global, immutable environment.
    ///
    /// This function does not check whether analytic update of the given time
    /// delta is actually possible. Generally, this is invoked with `et`
    /// increasing by one and then calling `step_nonanalytic()` alternately;
    /// however, the `Entity` may know in certain cases that the entity is pure
    /// analytic across some time-frame and skip this alternation, directly
    /// using the desired time delta.
    ///
    /// This step may change/advance the edges of the implied bounding box of
    /// this entity.
    fn step_analytic(&mut self, env: &Self::Env, et: Chronon);
    /// Performs a "non-analytic" step forward of this entity.
    ///
    /// A non-analytic step here is simply one that cannot be done in constant
    /// time for any delta or has an external effect (specifically, spawning
    /// new entities and/or ceasing to exist).
    ///
    /// Unlike analytic steps, nonanalytic steps are always done one chronon at
    /// a time.
    ///
    /// Note that, like `step_analytic()`, this function cannot observe other
    /// entities or other external mutable things, but simply the immutable
    /// environment.
    ///
    /// This step may *not* advance the edges of the implied bounding box of
    /// this entity, though it may relocate it if the bounding box shrinks to
    /// match.
    fn step_nonanalytic(&mut self, env: &Self::Env) ->
        EntityMutation<Self>;
    /// Updates this entity by applying the given event.
    ///
    /// This step may *not* advance the edges of the implied bounding box of
    /// this entity, though it may relocate it if the bounding box shrinks to
    /// match.
    fn apply_event(&mut self, event: &Self::Event, env: &Self::Env) ->
        EntityMutation<Self>;
    /// Returns the current position of this entity.
    fn position(&self) -> Position;
    /// Returns the size of the bounding box of this entity, measured from
    /// `position()` to each edge.
    fn size(&self) -> Dimension;
    /// Notifies this entity that it is an active participant in a collision
    /// with another entity `that`.
    ///
    /// Collisions do not necessarily take place at the entities' current
    /// position; if they are moving quickly enough relative to each other,
    /// collision may instead occur at a point between steps. `self_pos` and
    /// `that_pos` give the effective locations of the entities where the
    /// collision actually occurred.
    ///
    /// `points` is a non-empty slice of locations indicating the pairs of
    /// collision tree leaves which collided between the objects, where
    /// `target_mask` bits in `self`'s tree matched `subject_mask` bits in
    /// `that`'s tree. These coordinates are in the same intermediate space as
    /// `self_pos` and `that_pos`.
    ///
    /// This is a non-analytic update, and is the only way for multiple
    /// entities to interact with each other.
    ///
    /// This step may *not* advance the edges of the implied bounding box of
    /// this entity, though it may relocate it if the bounding box shrinks to
    /// match.
    fn collide(&mut self, self_pos: &Position,
               that: &Self, that_pos: &Position,
               points: &[CollisionLocation]) -> EntityMutation<Self>;
    /// Returns a reference to the root of the collision tree for this entity.
    fn collision_tree(&self) -> &CollisionTree;
    /// Discards any cached state held by this entity. This must have no
    /// observable effect on the entity's state.
    fn prune(&mut self);
}

/// Uniquely identifies an event within a single entity.
///
/// This is conceptually similar to the identifying event triple
/// (instant,object,serno) in Praefectus, but does not include the object id
/// since all events applied to one entity are necessarily part of the same
/// object.
///
/// Order is important --- events are applied to an entity in strictly
/// ascending order by the (time,serno) tuple.
#[derive(Copy,Clone,Debug,PartialOrd,Ord,PartialEq,Eq)]
pub struct EventId {
    /// The time at which an event applies.
    pub time: Chronon,
    /// Uniquefier to distinguish multiple events at the same instant.
    pub serno: praef::EventSerialNumber,
}

/// Identifies an entity within a praefectus system.
///
/// Ids are not stored on the entities themselves; they are implicit from
/// context. Their most important application is in dictating the total order
/// of collision application.
#[derive(Copy,Clone,Debug,PartialOrd,Ord,PartialEq,Eq)]
pub struct EntityId {
    /// The praefectus object which owns this entity.
    pub object: praef::ObjectId,
    /// The index of this entity within the object's entity array.
    pub index: i32,
}

/// Describes a subinstant of non-analytic entity transition.
///
/// This dictates the total order of non-analytic transitions within a single
/// instant: Non-analytic update, external events sorted ascending by serial
/// number, collisions sorted ascending by collidee, and spawning.
#[derive(Copy,Clone,Debug,PartialOrd,Ord,PartialEq,Eq)]
pub enum TransitionSub {
    /// A non-analytic local state change produced by
    /// `EntityType.step_nonanalytic()`.
    Update,
    /// The application of any external event.
    Event(praef::EventSerialNumber),
    /// A non-analytic state change produced by
    /// `EntityType.collide()`.
    Collision(EntityId),
    /// The transition of an entity from being non-existent to being existent.
    Spawn,
}

/// Uniquely identifies a point of non-analytic state transition within a
/// single entity.
///
/// This provides a total order for all non-analytic transitions.
#[derive(Copy,Clone,Debug,PartialOrd,Ord,PartialEq,Eq)]
pub struct TransitionId {
    /// The time of this transition.
    pub time: Chronon,
    /// The sub-instant of this transition, breaking ties with other
    /// transitions with the same time.
    pub sub: TransitionSub,
}

/// Trait for a value that can be stored in a tape.
///
/// It simply provides access to the underlying ordering key.
trait TapeValue: Eq {
    type Key: Ord + Eq;

    fn key(&self) -> &Self::Key;
}

/// Maintains a seekable sequence of "past" and "future" values.
///
/// There is at all times a "current" value, which is the value with the
/// greatest key which is less than or equal to the seek target.
#[derive(Clone)]
struct Tape<T: TapeValue> {
    content: Vec<T>,
    ix: usize
}

impl<T: TapeValue> Tape<T> {
    /// Creates a tape containing only the given "zero" element, which becomes
    /// the current value.
    fn new(zero: T) -> Tape<T> {
        Tape {
            content: vec![zero],
            ix: 0,
        }
    }

    /// Returns a reference to the current value in this tape.
    fn curr(&self) -> &T {
        &self.content[self.ix]
    }

    /// Returns the next value in this tape, or None if at the end.
    fn next(&self) -> Option<&T> {
        if self.ix + 1 < self.content.len() {
            Some(&self.content[self.ix + 1])
        } else {
            None
        }
    }

    /// Repositions the tape so that the "current" element is the one with the
    /// greatest key less than or equal to target.
    ///
    /// This is optimised for small movements, and uses a simple linear search.
    fn seek(&mut self, target: &T::Key) {
        while self.ix > 0 && *target < *self.curr().key() {
            self.ix -= 1;
        }
        while self.ix + 1 < self.content.len() &&
                *target >= *self.content[self.ix+1].key() {
            self.ix += 1;
        }
    }

    /// Adds a value to this tape.
    ///
    /// The tape must have been `seek()`ed to a key preceding the key of this
    /// value or otherwise be known to have that state.
    ///
    /// If the new value exactly matches the next future value, the new value
    /// is discarded, the old value becomes the new current, and true is
    /// returned, indicating that the future timeline is unchanged. Otherwise,
    /// any future values are discarded and the new value becomes the current,
    /// and false is returned to indicate that future information is now
    /// unknown.
    fn push(&mut self, value: T) -> bool {
        debug_assert!(*value.key() > *self.curr().key());

        if self.ix + 1 < self.content.len() &&
            value == self.content[self.ix + 1]
        {
            self.ix += 1;
            true
        } else {
            self.drop_future();
            self.push_direct(value);
            false
        }
    }

    fn push_direct(&mut self, value: T) {
        debug_assert!(self.ix + 1 == self.content.len());
        self.content.push(value);
        self.ix += 1;
    }

    /// Drops all future values for the tape.
    fn drop_future(&mut self) {
        self.content.truncate(self.ix + 1);
    }
}

const BEGINNING_OF_TIME: TransitionId = TransitionId {
    time: 0,
    sub: TransitionSub::Update,
};

#[derive(Clone, PartialEq, Eq)]
struct EntityVersion<T: Eq> {
    id: TransitionId,
    value: Option<T>,
}

impl<T: Eq> TapeValue for EntityVersion<T> {
    type Key = TransitionId;

    fn key(&self) -> &TransitionId {
        &self.id
    }
}

impl TapeValue for TransitionId {
    type Key = TransitionId;

    fn key(&self) -> &Self {
        self
    }
}

/// Manages the full state of an entity.
///
/// This is essentially a collection of known past and future non-analytic
/// transitions and a cached "current" state, plus all events applied to the
/// entity.
///
/// Entities at this level do not have dynamic lifetimes (though they may be
/// allocated on-demand); they notionally always exist, but may at times
/// represent non-existent entities, and undergo non-analytic transitions to
/// begin representing existent entities.
///
/// As hinted in the definition of `TransitionId`, the span of a single chronon
/// involves multiple discrete update steps. A chronon itself refers to the
/// state immediately before the analytic step. The sequence for each instant
/// is analytic, non-analytic, event, collision, spawn.
///
/// A fundamental state of an entity is whether it is currently "awake" or
/// "sleeping". A sleeping entity is one which currently knows its future
/// non-analytic transitions and has not been perturbed from the timeline which
/// leads to those transitions. Thus tests for non-analytic updates can be
/// completely elided; in particular, there is no need to test two sleeping
/// entities for collisions, because whether such a collision occurs is already
/// encoded in the future states. Wakeness is determined by the current `limit`
/// value of an entity relative to a desired time; if the `limit` is greater
/// than the sample time, the entity is sleeping at that time. The `limit`
/// recedes in response to retroactive changes, and is advanced by steping the
/// entity through its state machine.
///
/// The transition sequence and wakeness govern the interaction model with an
/// entity. When an entity is sleeping, it can simply be queried for its state
/// at any instant before its `limit`. Otherwise, the `limit` must be advanced
/// iteratively via the `step_internal()`, `collide()`, `spawn()`, and `next()`
/// functions.
///
/// If a sleeping entity encounters new information (for example, by colliding
/// with an awake object), this must be fed to it as with awake entities. This
/// may cause the entity to become awake, discarding future information and
/// receding `limit`, if the event perturbs its state to differ from what had
/// been calculated previously.
///
/// An `Entity` obeys the `SPEED_OF_LIGHT` restrictions except at certain
/// transitions called _warps_, which generally correspond to spawning events.
/// An `Entity` can be queried for the next instant at which it warps so that
/// this can be taken into account.
#[derive(Clone)]
pub struct Entity<T: EntityType> {
    /// The cached "current" state of this entity.
    ///
    /// This cannot be used as the analytic base; that must be retrieved from
    /// the end of `past`.
    current: EntityVersion<T>,
    /// Tape of past and future non-analytic transitions for this entity.
    transitions: Tape<EntityVersion<T>>,
    /// Tape of transitions at which warps occur.
    warps: Tape<TransitionId>,
    /// The first instant at which the non-analytic transitions for this entity
    /// are unknown.
    limit: Chronon,
    /// Map of events which are applied to this entity.
    events: BTreeMap<EventId,T::Event>,
}

impl<T: EntityType> Entity<T> {
    /// Creates a new Entity representing a non-existent entity.
    pub fn empty() -> Self {
        Entity {
            current: EntityVersion {
                id: BEGINNING_OF_TIME,
                value: None,
            },
            transitions: Tape::new(
                EntityVersion {
                    id: BEGINNING_OF_TIME,
                    value: None
                }),
            warps: Tape::new(BEGINNING_OF_TIME),
            limit: 0,
            events: BTreeMap::new(),
        }
    }

    /// Returns the time of the next warping transition, or None if this entity
    /// will not warp for the forseeable future.
    pub fn next_warp(&self) -> Option<TransitionId> {
        self.warps.next().map(|r| *r)
    }
}
