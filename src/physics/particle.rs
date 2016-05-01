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

/*!

Hexeline's physics model is somewhat unusual because it needs to be able to
efficiently step time forwards and backwards so as to apply events that arrived
late. It is further complicated in that it must be able to produce exactly the
same result regardless of the order in which this happens.

Terminology: A "particle" is a body with a two-dimensional position, a size,
and whose bounding box never advances at above the `SPEED_OF_LIGHT`. A particle
comes into existence and ceases to exist at definite points in time. A "pslot"
is a location which holds zero or one particle and manages the particle's
historical state. Each pslot belongs to one praefectus object and has a 32-bit
index which uniquely identifies it within that object.

Note that while a pslot manages a particle's historical state, we usually refer
to the particle itself more abstractly to include its state over time.

Let us start the design discussion by considering a system with only one
particle.

A simple design would be to advance the particle one step forward at a time,
saving each resulting state. Then, rewinding can simply jump directly to the
appropriate saved state. This is undesirable because it requires having coarse
time-steps or requires a lot of memory for the intermediate states.

The first optimisation is that, in many cases, the future state of a particle
for many steps forward can be determined from a prior state for the same cost
as advancing one step. For example, updating position based on velocity. Thus,
we split updates into "analytic" updates where this is possible, and
"non-analytic" updates where it is not. Because we can do any number of
analytic steps forward at once, we only need to save states after non-analytic
state transitions. Furthermore, if we know when the next non-analytic state
transition is, we don't even need to check for non-analytic updates.

Note, however, that since analytic updates over N steps may be an approximation
of doing N individual steps, we must *only* make snapshots at non-analytic
state transitions. Any state after a non-analytic transition is derived by
doing one analytic update for the appropriate time differential from that
non-analytic transition, rather than continuing from some "current" state
derived from a prior analytic update.

Generally, when a historical version of a particle is recalled, we still know
the future non-analytic transitions, and can thus jump cheaply forward as well
up to some limit, the "deterministic limit". In a one-particle system, the
limit can be recessed by the addition or removal of external events (from
praefectus). When this happens, the deterministic limit is reset to the instant
of the affected event and all future states beyond the deterministic limit are
dropped. At a given instant T, a particle whose deterministic limit is greater
than T is "passive"; a particle whose deterministic limit is less than or equal
to T is "active".

Things become much more difficult once the system has multiple particles that
can interact with each other. Notionally, such interactions occur via
observations and via collisions.

Observations are particularly tricky, as a particle could make seemily
arbitrary transitions in response to the exact state of any other particle in
the system. This would require all particles to be advanced in lock-step
whenever *any* particle needs to be tested for non-analytic transitions,
largely defeating the purpose of all we've done. Instead, we simply ban direct
observation; this concept is instead represented via collision, so that any
solution for collision equally applies to observation.

The concept of "collision" needs to be broadened as a result.

- There are multiple "types" of collision. For example, physical and
  observational collision are different. Particles may have different collision
  boundaries for different collision types.

- Collision is not reflexive. Ie, if A collides with B, B does not necessarily
  collide with A. More strictly, if B does not collide with A, regardless of
  whether A collides with B, B's state is not directly affected by A's state at
  that time. In tfhe statement "A collides with B", A is referred to as the
  "observer" and B as the "subject" ("collider"/"collidee" are not used since
  they are easily confused).

The statement "A collides with B at instant T" is specifically defined to mean
that A may change its state based on B's exact state at instant T. However,
collision state changes must be commutative and associative. This allows the
collision detection system to prioritise efficiency.

Each particle defines a tree of axis-aligned bounding boxes which define its
collision boundaries. Each bounding box has an observer mask and a subject
mask; one box collides with another only if they are in different particles and
the intersection of the observer's observer mask and the subject's subject mask
is non-empty.

No edge of any bounding box of any particle may advance at a speed in exceess
of the `SPEED_OF_LIGHT`. Note that there is no such restriction for the actual
position of a particle provided its bounding boxes reshape in a way that does
not violate the speed limit, nor is there any limit on the rate at which a
bounding box may recess. The observer/subject masks on collision boxes
generally cannot be augmented post-facto, since this would have the effect of a
collision boundary instantly expanding from nothing.

These properties mostly get us to where we want:

- Passive particles capture the results of collisions within their future
  states. This collisions between passive particles need not be recalculated;
  only collisions involving at least one active particle must be considered.

- The exact state of a passive particle at a certain instant is often not
  required to test for collision with an active particle, as the
  `SPEED_OF_LIGHT` limit allows extrapolating whether collision is even
  theoretically possible.

What we still need is propagation of nondeterminism. One aspect of this is
obvious: If passive particle A collides with active particle B at instant T,
A's deterministic limit is reduced to T, making A active.

A more subtle issue is that if passive particle A's future state captures a
collision with active particle B at instant T, but A does not in fact collide
with B at that instant, A's deterministic limit must also be reduced to T. This
we accomplish by tracking the greatest instant at which each pslot has been the
subject of collision with each other pslot. When pslot A's deterministic limit
reduces to T, all pslots whose last-collision time in A's table is greater than
or equal to T have their deterministic limit reduced to T as well. (This is
more conservative than it could be, but the smaller benefit of more precisely
tracking the exact instants and removing non-collisions would likely be
counter-balanced by the cost of extra bookkeeping.)

The last awkward thing to deal with is the creation and destruction of
particles, which are necessarily mapped into shared pslots. The first question
is whether pslot assignment needs to be deterministic --- ie, whether all nodes
in the system must place the same particles into the same pslots. Since
collisions are unordered, we do not need particles to be ordered either, so as
long as pslot ids are not sent over the network or used for state transitions,
they do not need to be deterministic.

It might initially seem desirable to be able to address pslots over the
network, to describe events directed at particular particles. This is not
actually a workable idea, because even deterministic pslots wouldn't be
fully deterministic, so events could be sent to the incorrect particle.
Instead, we simply assign ids to particles that may be interested in events,
which are provided explicitly in the events that spawn the interested
particles.

There are three causes of existence changes:

- An external event (not applied to a particular particle) spawns a particle.

- A non-analytic update to a particle spawns another particle.

- A partle ceases to exist as a result of its non-analytic update.

The third case by itself is fairly simple, as it is simply a non-analytic state
change for the pslot as a whole. The awkwardness comes in with how it interacts
with the other two.

The first case can be merged into the second case by having a virtual particle
which always exists and is the target of such events, and thus spawns those
particles in its own non-analytic transition. This virtual particle is always
in pslot 0 and is called the "mother particle".

The spawning of one particle by another can be modelled as a mutual collision
for bookkeeping purposes. The parent pslot observes one or more other pslots to
find an empty one; thus, if any of those pslots become active before the spawn
instant, the parent needs to become active again since it may need to assign
different pslots. The child observes the parent since changes to the parent's
history could affect whether the child is actually spawned or which pslots are
used, etc.

The collision detector needs to be aware of the next spawn event for each
pslot, because spawning a particle effectively violates the `SPEED_OF_LIGHT`
limit. This is represented as a "warp", simply indicating an instant at which a
particle violates the speed of light. We also allow normal nonanalytic updates
to perform warps, as it simplifies a number of things such as adding cells to
ships.

The order for updates to a pslot/particle in each instant is as follows:

- Analytic update from the last non-analytic transition.

- (Active only) Apply collisions in any order. Instant is non-analytic if state
  changes.

- (Active only) Apply external events in order. Instant is non-analytic if
  state changes.

- (Active only) Check for non-analytic update. Particle spawning and
  nonexistence may happen here. Instant is non-analytic if state changes,
  particles are spawned, or particle ceases to exist.

The state after all steps is saved if the instant for that pslot is determined
to be non-analytic.

Collision application must be the first potentially non-analytic step as it is
the one step that can cause a particle to become active during the update
proper. Were that step elsewhere, we'd need to store intermediate non-analytic
states as well to properly step to the appropriate state. (Deterministic pslot
spawning can cause another pslot's deterministic limit to recess, but not to
the instant being processed.)

 */

#![allow(dead_code)]

use std::fmt::Debug;

use smallvec::SmallVec;

use super::defs::*;

/// Type for the collision observer/subject masks.
pub type CollisionMask = u16;
/// Type for indexing pslots
pub type PslotIndex = u32;

/// A node in a bounding box tree.
///
/// A node may be a leaf node or a branch node. Branch nodes themselves do not
/// result in detected collisions; rather, the branch is split into its
/// children and each child is tested for collision instead.
pub struct CollisionTreeNode<'a> {
    /// The mask of collision types for which this node can be an observer.
    pub observer_mask: CollisionMask,
    /// The mask of collision types for which this node can be a subject.
    pub subject_mask: CollisionMask,
    /// The position of the centre of this bounding box.
    pub pos: Position,
    /// The X and Y dimensions of the bounding box, measured from the centre to
    /// the edge (so they measure something akin to "radius" in L-infinity
    /// space).
    pub radius: Dimension,
    /// If present, this is a branch node, and the given thunk can be invoked
    /// to obtain the node's children. Otherwise, this node is a leaf node.
    pub children: Option<&'a CollisionTreeThunk>,
    /// Integer associated with collision results. Only useful on leaf nodes.
    /// The exact interpretation depends on the particle to which this node
    /// belongs.
    pub id: u32,
}

/// Thunk used to obtain the children for a collision tree node.
///
/// This exists so that particles don't need to generate the full tree every
/// frame, since in most cases only the top nodes are used, if any.
pub trait CollisionTreeThunk {
    fn get(&self) -> &[CollisionTreeNode];
}

/// Represents a pair of points in two particles which collided.
#[derive(Copy,Clone,Debug)]
pub struct CollisionPoint {
    /// The id on the node of the observer
    pub observer_id: u32,
    /// The id on the node of the subject
    pub subject_id: u32,
    /// The position of the observer collision tree node
    pub observer_pos: Position,
    /// The position of the subject collision tree node
    pub subject_pos: Position,
    /// The bitwise AND of the observer's observer mask and the subject's
    /// subject mask. Always non-zero.
    pub types: CollisionMask,
}

/// The result of a nonanalytic step which resulted in a nonanalytic state
/// transition.
#[derive(Clone,Debug)]
pub struct NonanalyticTransition<T: Clone + Debug> {
    /// Whether the particle that just updated still exists.
    pub exists: bool,
    /// Whether this particle "warped", ie, it may have expanded a collision
    /// boundary faster than the speed of light, or created new collision boxes
    /// in the tree that may have had the same effect.
    pub warp: bool,
    /// Other particles spawned as a result of this transition.
    pub spawn: SmallVec<[T;1]>,
}

/// Trait describing a particle. See the module documentation for details.
pub trait Particle: Clone + Debug {
    /// Type for the static, immutable "environment" data all particles can
    /// observe freely.
    type Env;
    /// Type for external events applied to particles.
    type Event;
    /// Type for holding heavier cached data about a particle which can be
    /// regenerated when needed and which should not be stored with historic
    /// versions of the particle.
    ///
    /// No guarantees are made about the identity, durability, or consistency
    /// of the cache. Particles must be prepared for the cache contents to be
    /// lost or rolled back at any time, and for caches belonging to other
    /// particles to be provided instead.
    type Cache : Default;

    /// Returns the calculated root of the collision tree for this particle.
    fn collision_bounds<'a>(&'a self, cache: &'a Self::Cache)
                            -> &'a CollisionTreeNode;
    /// Returns whether the given event is "addressed to" this particle.
    ///
    /// This must be unaffected by calls to `advance_analytic()`.
    fn is_addressee_of(&self, event: &Self::Event) -> bool;

    /// Performs an analytic advance of the particle by the given number of
    /// steps.
    fn advance_analytic(&mut self, cache: &Self::Cache,
                        env: &Self::Env, delta: Chronon);
    /// Returns a lower bound number of steps this particle can be advanced via
    /// advance_analytic() before a nonanalytic transition will be necessary,
    /// assuming no interactions with external events or other particles occur.
    /// Must return at least 1.
    fn max_analytic_advance(&self, cache: &Self::Cache, env: &Self::Env) -> u16;
    /// Notifies the particle that it collided with another particle at the
    /// given collision points.
    ///
    /// This function must be commutative and associative, and must not depend
    /// on the order of the points slice. Note that this also applies to
    /// collisions applied to the `other` particle; effectively, this means
    /// that this call cannot depend on data that `other.collide_with()` might
    /// itself mutate.
    ///
    /// Returns whether the particle underwent any state change. If this
    /// returns false, the collision is treated as if it did not occur.
    fn collide_with(&mut self, self_cache: &mut Self::Cache,
                    other: &Self, other_cache: &Self::Cache,
                    env: &Self::Env,
                    points: &[CollisionPoint]) -> bool;
    /// Applies the given event to this particle. Only events which passed
    /// `is_addressee_of()` are passed to this function.
    fn apply_event(&mut self, cache: &Self::Cache,
                   event: &Self::Event, env: &Self::Env);
    /// Performs a non-analytic step on this particle.
    ///
    /// If the particle changes state in any way, returns a structure
    /// describing that state change and any particles spawned as a result.
    ///
    /// If the particle does not change state at all, returns None.
    fn advance_nonanalytic(&mut self, &mut Self::Cache)
                           -> Option<NonanalyticTransition<Self>>;
}
