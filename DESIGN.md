HIGH-LEVEL GAMEPLAY

Initial player configuration: Selected race.

A player starts out with a small but functional basic ship for their race
(randomly chosen from a set of presets) and a certain amount of free materials.
Game modes initially give the player some form of peaceful time to spend these
materials if desired (be this through invincibility or lack of enemies).

Players are free to fly about the map and attack whatever they please.
Destroying things generally releases materials the player may pick up. At any
time, a player can spend materials to add components to their ship.
Additionally, they may voluntarily destroy a component on their ship to get
half its constituent materials. These actions are usually done in realtime, but
easier single-player modes may pause instead.

Damage is taken on a component-by-component basis. Different components are
sensitive to different types of damage, and may damage other things nearby when
destroyed. If a ship has taken no damage for a set amount of time, all
components heal damage. Destroyed components do not get restored.

BASIC CONTROLS

WASD / Arrow keys: Accelerate in given directions. "Backwards" is usually
"brake". Most engine types can only accelerate in one direction.

Shift: Overdrive engines

Mouse left/right: Rotate ship.

Left click: Fire current weapon

Right click: Tractor nearby free materials / other effects

Scroll: Change weapons (as in Sauerbraten, but with a display on the screen so
newcomers can understand it)

Space: Change target

C/V: Cycle camera

Tab: Open ship management console

RACES

Terrans:

- No upsides

- No downsides.

White/bright red/yellow.

The official military and paramilitary forces of unified Earth.

Pirates:

- Higher-level components are strictly superior to lower-level components.

- Only base-level components immediately available. Higher levels acquired via
  upgrades after damage/repair cycles on existing components.

Dark grey/mute violet/white

The remnants of a rebellion against the unified Earth nation which now resides
in a semi-organised state in the Centauri system. Their technology is decades
out of date, but they remain a formidable force through sheer tenacity.

Replicants:

- No native higher-level components.

- Can construct any component on any live visible ship.

Green/black/yellow.

A mostly non-intelligent family of self-replicating military robots introduced
by one of the losing contenders for governance of unified Earth. Their own
technology is half a century obsolete, but they have the unique ability to
replicate new technology they see.

Marians:

- -50% hull damage; -75% shield damage; +25% weapon damage; -75% shock damage.

- -75% weapon fire rate (weapon damage adjusted to cancel out). -50% engine
  thrust. Dynamic shields unavailable and spread weapons unavailable.

Blue/dark blue/dark cyan.

An ancient aquatic race not known for doing anything particularly quickly.
Their ships are built to endure and achieve victory via attrition rather than
surprise attack. Shock damage is relatively ineffectual since it mostly simply
counter-balances the immense liquid pressure the hull must wistand.

Jovians:

- +100% weapon damage, -50% materials costs. +100% engine thrust.

- -50% hull strength, +100% heat generation. Energy shields unavailable.

Yellow/black/red

A reptilian race of unknown origin which was first discovered after having
established a base in orbit of Jupiter, which set into motion the events that
eventually led to Earth's unification. They generally value cheap construction
and high destructive power over safety.

Striders:

- Equally good acceleration in all directions. All weapons guided.

- Momentum not conserved; ship quickly deccelerates when engines not running.
  Weapon guidance must be assisted with energy pulses emitted manually (ie,
  with secondary fire). No heavy weapons.

Light grey/black/white

Little is known about this mysterious species, except that their engines and
guidance systems appear to work by dragging space time itself.

MATERIALS

Constructing new ship components consumes resources; each component type has
different resource requirements. Every resource other than Tech can be found
naturally in space. Additionally, destroyed components release 50% of their
constituent materials into space.

Each ship has a fixed maximum capacity for unused materials; this forces the
player to invest in new construction early and incrementally improve their
ship, rather than stockpiling resources for a sudden massive change to the
ship.

When a ship is neutralised, all unused materials it is carrying are released
into space.

Alloy: Metals for constructing hull, system enclosures, and so forth. The most
common resource. Obtainable from astroids and general debris.

Conductor: More refined metals and ceramics which are good conductors of
eletricity. Occasionally found in asteroids.

Tech: Materials in general which do not occur naturally. The bridge and many
higher-level component types contain Tech.

SHIP/COMPONENT SYSTEM

A ship is composed of a number of components, all of which are equally-sized
hexagons. All components in a ship are connected; some components have
limitations regarding how they can be connected and whether certain paths need
to be infinitely unobstructed.

Components are grouped into the following categories:

- Bridge. Every ship has exactly one bridge component.

- Weapons (subdivided into direct/spread/missile/melee)

- Energy sources (subdivided into electric/anti-heat/hydrogen)

- Energy storage (subdivided as with energy sources)

- Energy transport (weakly subdivided as with energy sources)

- Shields and armour

- Engines and boosters

Energy is a key component in the functioning of a ship. Most components require
at least one supply of energy to function. Energy is produced at energy
sources. Generally, energy is only delivered to the immediate neighbours of an
energy source. Energy transport components are used to extend the range. A
notable exception is anti-heat, which is conducted at a low rate across all
components and slowly added to all components over time. There are three energy
types:

- Electric. The primary energy type which powers almost all components.
  Components which do not receive sufficient electric energy simply do not
  function in most cases, but certain things may fail catastrophically (eg, a
  plasma burst launcher interrupted mid-charge would release its accumulated
  plasma directly upon itself and surrounding components).

- Anti-heat. Many components consume anti-heat (ie, generate heat) during their
  functioning. All components store anti-heat and conduct it to at least some
  degree. A component will be destroyed (possibly violently) if its anti-heat
  becomes too low. Some components function differently at different levels of
  stored anti-heat. "Heat energy" is negated so that the numbers here remain
  positive for all energy types.

- Hydrogen. Atomic hydrogen gas collected from the stellar medium. Required by
  certain engine types and weapons, and can be used to augment certain
  components' functionality.

Energy transport components are unique in that they may be stacked up to 3
deep, so that wire-crossing is not an issue. Each energy transport component
has a maximum flow rate it can sustain.

The requirement to be able to *transport* energy from source to sink is
designed to counteract a central scaling issue seen prominently in both
Abendstern and Reassembly: Increasing the ship's radius linearly only increases
attack surface linearly, but increases power supply and damage output
quadratically. In this system, a player has a spectrum of choices between two
extremes:

- Keep energy sources in a safe place and add infrastructure to move the
  energy. Increasing energy supply quadratically also requires quadratically
  increasing the *width* of any energy transport. However, the length of each
  transport is linearly proportional to the radius of the ship, so the total
  area consumed by transport infrastructure grows cubically. Thus, centralised
  designs cannot be scaled to very large ships.

- Make the ship homogenous: energy sources must be colocated with their sinks.
  First of all, this makes the ship much more vulnerable, as components which
  fail violently must be near the surface. Secondly, the ship must either have
  a bunch of isolated systems with low peak power, or be a large interconnected
  system which allows cascading failures to spread quickly.

The standard electrical terms of energy (joules), charge (coulombs), potential
(voltage), current (amperage), power (wattage), resistance (ohms), and
capacitance (farads) are used to describe all three energy types of energies.
They are thus qualified, eg, electric-volts (eV), anti-heat-ohms (hΩ), or
hydrogen-watts (HW).
