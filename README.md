# Kepler Orbit - 3D Apsidal & Nodal Precession Demo

An interactive 3D simulation (VPython) of the restricted two-body
(Kepler) problem, extended with two real orbital-precession effects:
**apsidal precession** (the periapsis itself rotates, like Mercury's
famous perihelion advance) and **nodal precession** (the whole orbital
plane slowly sweeps around the polar axis, like the regression caused
by a planet's oblateness). You choose the astronomical parameters --
including how strong you want each precession effect to be -- in the
terminal before the 3D scene opens.

![orbital precession](https://upload.wikimedia.org/wikipedia/commons/2/2f/Relativistic_precession.svg)

---

## Why this is a genuinely good physics demo

This isn't a cosmetic "wobbly ellipse" effect bolted on for visual
interest -- both kinds of precession are implemented as the actual
physics that produce them in nature, deliberately built so you can see
the underlying mechanism, not just its result:

- **Apsidal precession from a real post-Newtonian correction.** The
  radial acceleration used in the simulation is

  ```
  a_r(r) = -GM / r^2  -  3 * GM * L^2 / (c^2 * r^4)
  ```

  where `L` is the (conserved) specific angular momentum. This second
  term is the standard weak-field General Relativity correction to
  Newtonian gravity -- the same one Einstein used to explain the
  anomalous advance of Mercury's perihelion, one of the first
  confirmed predictions of GR. It predicts a precession per orbit of

  ```
  delta_phi = 6 * pi * GM / (c^2 * a * (1 - e^2))
  ```

  The script inverts this formula: you say how many degrees of
  precession you want per orbit, and it solves for the effective `c^2`
  needed to produce that. While the simulation runs, it *independently
  measures* the precession the numerical integration actually
  produces (by tracking successive periapsis passages) and displays it
  live next to your target. Because the formula above is only the
  leading-order approximation, the two numbers agree closely for
  modest, realistic precession values and increasingly diverge for
  large, deliberately exaggerated ones -- a nice hands-on illustration
  of first-order vs. full nonlinear behavior.

- **A symplectic integrator, chosen for a physical reason.** The
  radial equation of motion is integrated with velocity Verlet rather
  than, say, RK4. Symplectic integrators conserve angular momentum and
  (on average) energy far better over long integrations. That matters
  a lot here specifically because the whole point is to measure a real
  physical precession signal -- an integrator that leaks energy or
  angular momentum would introduce its *own* spurious numerical
  precession and contaminate the measurement. (You can check this
  yourself: set the desired precession to 0 and watch the "measured"
  readout -- with a correctly conserving integrator it stays at
  essentially zero, confirming the ellipse doesn't drift on its own.)

- **Nodal precession as a genuinely three-dimensional effect.** Real
  orbits also precess out of their own plane -- a planet's equatorial
  bulge (the J2 effect) or, far more subtly, relativistic frame
  dragging (Lense-Thirring precession) causes the ascending node to
  regress over time. The simulation applies this as a user-chosen rate
  in degrees per orbit, rotating the whole orbital plane about the
  polar axis while the in-plane rosette from apsidal precession keeps
  evolving inside it -- which is exactly what makes the visualization
  actually three-dimensional rather than a flat ellipse viewed from an
  angle.
- **Boltzmann-free, but the same "physics-first" spirit.** Like a
  simulated-annealing optimizer borrows the Boltzmann factor from
  statistical mechanics, this demo borrows the exact perturbative
  correction General Relativity makes to Newtonian gravity -- a good
  companion piece for anyone using these scripts to show how directly
  physical laws translate into working, visualizable algorithms.

---

## Installation and running

```bash
git clone https://github.com/<your-username>/kepler-precession-3d.git
cd kepler-precession-3d

python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

pip install -r requirements.txt

python kepler_precession_3d.py
```

You will be prompted in the terminal for:

1. **GM** -- the gravitational parameter of the central body (`G *
   mass`, combined into one number, as is standard practice in
   astrodynamics since the orbital dynamics only ever depend on this
   product).
2. **Periapsis distance** -- the closest approach of the orbit.
3. **Eccentricity** -- `0` for a circle, up to `0.9` for a strongly
   elongated ellipse.
4. **Desired apsidal precession per orbit** (degrees) -- `0` for a
   plain, non-precessing Newtonian ellipse.
5. **Orbital inclination** (degrees) -- `0` for flat in the reference
   plane, `90` for a polar orbit.
6. **Nodal precession rate** (degrees per orbit) -- can be negative;
   `0` disables it.
7. **Number of orbits to simulate.**

Press Enter at any prompt to accept the default shown in brackets.

Once the 3D scene opens: drag with the mouse to rotate the view,
scroll to zoom. Two buttons below the scene let you **pause/resume**
the animation and **hide/show** the trajectory trail. A live status
line reports the current orbit count, the current node orientation,
and the target vs. measured apsidal precession.

The **gray ring** marks the reference (equatorial) plane, the **faint
static curve** is the unperturbed, non-precessing Kepler ellipse for
comparison, and the **bright curve** is the actual, evolving simulated
trajectory.

---

## How it works

- **State:** polar coordinates `(r, theta)` within the orbital plane.
  Using polar coordinates keeps `theta` continuous and unwrapped
  (no branch-cut headaches) and makes periapsis detection --
  local minima of `r(t)` -- straightforward.
- **Radial equation:** `r'' = a_r(r) + L^2 / r^3`, integrated with
  velocity Verlet (see above for why).
- **Angular equation:** `theta' = L / r^2` (Kepler's second law /
  conservation of areal velocity), integrated with the trapezoidal
  rule.
- **Initial conditions:** the satellite starts at periapsis, where the
  velocity is purely tangential; periapsis speed follows directly from
  the vis-viva equation given the periapsis distance and eccentricity,
  so there is no need to specify a velocity vector by hand.
- **3D embedding:** the in-plane `(r, theta)` solution is rotated into
  3D using the standard perifocal-to-inertial transform for a fixed
  inclination and a continuously growing longitude of the ascending
  node (the nodal precession).
- **Precession measurement:** successive periapsis passages are
  located with a parabolic interpolation of `r(theta)` near each local
  minimum, and the drift in periapsis angle between consecutive
  passages is averaged into the "measured precession per orbit"
  readout.

Fixed constants such as the animation frame rate, target animation
length, and integration resolution can be tuned at the top of
`kepler_precession_3d.py`.

---

## Repository structure

```
kepler-precession-3d/
├── kepler_precession_3d.py   # the whole simulation: physics + 3D visualization
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Possible improvements

- Add a second, non-precessing "ghost" satellite running in parallel
  (analogous to the shadow trajectory in a chaos demo), so the
  precessing and non-precessing orbits can be compared directly, side
  by side, in real time.
- Derive the nodal precession rate from an actual J2 oblateness
  parameter of the central body instead of setting it directly, for a
  more physically grounded (if less immediately controllable) version
  of the effect.
- Plot the measured precession per orbit over time to show how quickly
  it converges to a stable average.

---

## License

Released under the MIT License (https://github.com/Anatase666/Kepler-orbit-simulation?tab=MIT-1-ov-file#) — see the file for details.
