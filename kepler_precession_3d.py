"""
Interactive 3D visualization of the Kepler (restricted two-body) problem,
extended with two kinds of orbital precession that are staples of real
astronomy and astrophysics:

1. Apsidal precession (rotation of the periapsis within the orbital
   plane), implemented with the standard first-order post-Newtonian
   correction to Newtonian gravity:

       a_r(r) = -GM / r^2  -  3 * GM * L^2 / (c^2 * r^4)

   where L is the (conserved) specific angular momentum. This is the
   same weak-field correction that famously explains the anomalous
   precession of Mercury's perihelion in General Relativity, and it
   predicts a precession per orbit of

       delta_phi = 6 * pi * GM / (c^2 * a * (1 - e^2))

   This script inverts that formula: you choose how many degrees of
   precession you want per orbit, and the script solves for the
   effective "c^2" needed to produce it -- then, as the simulation
   runs, it independently *measures* the precession actually produced
   by the numerical integration (by tracking successive periapsis
   passages) and reports it next to the target. Since the analytic
   formula above is only the leading-order (small-precession)
   approximation, the two numbers agree closely for modest precession
   values and increasingly diverge for large, deliberately exaggerated
   ones -- itself a small illustration of first-order vs. full
   nonlinear behavior.

2. Nodal precession (rotation of the whole orbital plane about the
   polar axis), of the kind real orbits experience from a planet's
   oblateness (the J2 effect) or, far more weakly, from relativistic
   frame dragging (Lense-Thirring precession). Here it is applied as a
   user-specified rate in degrees per orbit, which is what makes the
   visualization genuinely three-dimensional: an inclined orbit's
   whole plane slowly sweeps around the polar axis while the orbit
   itself traces out a precessing rosette within that plane.

Physics and integration:
    The motion is integrated in polar coordinates (r, theta) in the
    orbital plane, which keeps the precession angle theta continuous
    and unwrapped (no branch-cut headaches) and makes periapsis
    detection straightforward. The radial equation

        r'' = a_r(r) + L^2 / r^3

    is integrated with velocity Verlet (a symplectic integrator),
    deliberately chosen over e.g. RK4: symplectic integrators conserve
    angular momentum and (on average) energy far better over many
    orbits, which matters a lot here since we are trying to measure a
    real, physical precession signal and do not want the integrator
    itself to introduce spurious numerical precession that would
    contaminate it.

    The angle theta is advanced with the standard areal-velocity
    relation theta' = L / r^2 (Kepler's second law), using trapezoidal
    integration.

    The 2D (r, theta) solution -- which already includes whatever
    apsidal precession the physics produces -- is then embedded in 3D
    using the standard perifocal-to-inertial rotation for a fixed
    inclination i and a continuously growing longitude of the
    ascending node Omega(t) (the nodal precession).

Controls:
    Parameters (gravitational parameter, periapsis distance,
    eccentricity, desired apsidal precession, inclination, nodal
    precession rate, number of orbits) are entered interactively in
    the terminal before the 3D scene opens; press Enter at any prompt
    to accept the default shown in brackets.

    Once the simulation is running, two buttons below the 3D scene let
    you pause/resume the animation and hide/show the trajectory trail.
    The scene itself is a normal VPython 3D view: drag with the mouse
    to rotate, scroll to zoom.
"""

import math

from vpython import (
    button,
    canvas,
    color,
    curve,
    cylinder,
    ring,
    sphere,
    vector,
    rate,
    wtext,
)

# --- Animation pacing (edit here for a faster/slower playback) ----------
RENDER_FPS = 60
TARGET_ANIMATION_SECONDS = 30.0
N_STEPS_PER_ORBIT = 720
MAX_TRAIL_POINTS = 20000

# --- Default astronomical parameters (scene units) -----------------------
DEFAULT_GM = 6.0e4
DEFAULT_R_PERIAPSIS = 120.0
DEFAULT_ECCENTRICITY = 0.6
DEFAULT_PRECESSION_DEG = 12.0
DEFAULT_INCLINATION_DEG = 25.0
DEFAULT_NODAL_RATE_DEG = 8.0
DEFAULT_NUM_ORBITS = 10


# --- Interactive parameter entry ------------------------------------------

def prompt_float(message, default, min_value=None, max_value=None):
    """Prompts for a float; an empty input falls back to `default`."""
    while True:
        raw = input(f"{message} [{default}]: ").strip()
        if not raw:
            return default
        try:
            value = float(raw)
        except ValueError:
            print("Please enter a valid number.")
            continue
        if min_value is not None and value < min_value:
            print(f"Please enter a number >= {min_value}.")
            continue
        if max_value is not None and value > max_value:
            print(f"Please enter a number <= {max_value}.")
            continue
        return value


def prompt_int(message, default, min_value=None, max_value=None):
    """Prompts for an integer; an empty input falls back to `default`."""
    while True:
        raw = input(f"{message} [{default}]: ").strip()
        if not raw:
            return default
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a valid whole number.")
            continue
        if min_value is not None and value < min_value:
            print(f"Please enter a number >= {min_value}.")
            continue
        if max_value is not None and value > max_value:
            print(f"Please enter a number <= {max_value}.")
            continue
        return value


def gather_parameters():
    print("=== Kepler Orbit -- 3D Precession Demo ===")
    print("(press Enter to accept the default value shown in brackets)\n")

    GM = prompt_float(
        "Gravitational parameter GM of the central body (G * mass)",
        DEFAULT_GM,
        min_value=1e-6,
    )
    r_periapsis = prompt_float(
        "Periapsis distance (closest approach)", DEFAULT_R_PERIAPSIS, min_value=1e-6
    )
    eccentricity = prompt_float(
        "Orbital eccentricity (0 = circle, closer to 1 = more elongated)",
        DEFAULT_ECCENTRICITY,
        min_value=0.0,
        max_value=0.9,
    )
    precession_deg = prompt_float(
        "Desired apsidal precession per orbit, degrees (0 = pure Newtonian ellipse, no precession)",
        DEFAULT_PRECESSION_DEG,
        min_value=0.0,
        max_value=90.0,
    )
    inclination_deg = prompt_float(
        "Orbital inclination, degrees (0 = flat in the reference plane, 90 = polar)",
        DEFAULT_INCLINATION_DEG,
        min_value=0.0,
        max_value=180.0,
    )
    nodal_rate_deg = prompt_float(
        "Nodal precession rate, degrees per orbit (rotation of the whole orbital "
        "plane, like J2 oblateness or frame dragging; 0 = none, may be negative)",
        DEFAULT_NODAL_RATE_DEG,
        min_value=-60.0,
        max_value=60.0,
    )
    num_orbits = prompt_int(
        "Number of orbits to simulate", DEFAULT_NUM_ORBITS, min_value=1, max_value=200
    )

    return {
        "GM": GM,
        "r_periapsis": r_periapsis,
        "eccentricity": eccentricity,
        "precession_deg": precession_deg,
        "inclination_deg": inclination_deg,
        "nodal_rate_deg": nodal_rate_deg,
        "num_orbits": num_orbits,
    }


# --- Orbital mechanics ------------------------------------------------------

def semi_major_axis(r_periapsis, eccentricity):
    return r_periapsis / (1 - eccentricity)

def periapsis_speed(GM, r_periapsis, eccentricity):
    """Speed at periapsis, purely tangential (Newtonian vis-viva relation)."""
    return math.sqrt(GM * (1 + eccentricity) / r_periapsis)

def orbital_period(GM, a):
    return 2 * math.pi * math.sqrt(a ** 3 / GM)

def c_squared_for_target_precession(GM, a, eccentricity, precession_rad):
    """
    Inverts delta_phi = 6*pi*GM / (c^2 * a * (1-e^2)) to find the effective
    c^2 that should produce the requested apsidal precession per orbit.
    Returns None for zero precession (i.e. "skip the correction term").
    """
    if precession_rad <= 0:
        return None
    return 6 * math.pi * GM / (precession_rad * a * (1 - eccentricity ** 2))

def radial_double_dot(r, GM, L, c_squared):
    """r'' = physical radial acceleration + centrifugal term (polar coords)."""
    a_r = -GM / r ** 2
    if c_squared is not None:
        a_r += -3 * GM * L ** 2 / (c_squared * r ** 4)
    return a_r + L ** 2 / r ** 3

def step_state(r, r_dot, theta, dt, GM, L, c_squared):
    """One velocity-Verlet step for r, r_dot, plus trapezoidal theta update."""
    a_now = radial_double_dot(r, GM, L, c_squared)
    r_next = r + r_dot * dt + 0.5 * a_now * dt * dt
    a_next = radial_double_dot(r_next, GM, L, c_squared)
    r_dot_next = r_dot + 0.5 * (a_now + a_next) * dt
    theta_next = theta + 0.5 * (L / r ** 2 + L / r_next ** 2) * dt
    return r_next, r_dot_next, theta_next


# --- 3D embedding: perifocal (in-plane) coordinates -> inertial frame ------

def perifocal_to_inertial(x_pf, y_pf, inclination, raan):
    """
    Standard perifocal -> inertial rotation (argument of periapsis omega = 0,
    since periapsis rotation is already captured by theta in the plane).
    """
    cos_o, sin_o = math.cos(raan), math.sin(raan)
    cos_i, sin_i = math.cos(inclination), math.sin(inclination)
    x = x_pf * cos_o - y_pf * cos_i * sin_o
    y = x_pf * sin_o + y_pf * cos_i * cos_o
    z = y_pf * sin_i
    return vector(x, y, z)


def reference_ellipse_curve(a, eccentricity, inclination, raan0, segments=200):
    """
    A static, faint curve showing the unperturbed (non-precessing) Kepler
    ellipse for comparison against the actual, precessing trajectory.
    """
    points = []
    for i in range(segments + 1):
        theta = 2 * math.pi * i / segments
        r = a * (1 - eccentricity ** 2) / (1 + eccentricity * math.cos(theta))
        x_pf, y_pf = r * math.cos(theta), r * math.sin(theta)
        points.append(perifocal_to_inertial(x_pf, y_pf, inclination, raan0))
    return points


# --- Scene setup -------------------------------------------------------------

def setup_scene(apoapsis):
    scene = canvas(
        title="Kepler orbit -- apsidal & nodal precession (drag to rotate, scroll to zoom)",
        width=1000,
        height=750,
        background=color.black,
        up=vector(0, 0, 1),
        forward=vector(-1, -0.6, -0.4),
        range=apoapsis * 1.3,
    )
    return scene


def draw_reference_geometry(apoapsis, central_radius):
    central_body = sphere(
        pos=vector(0, 0, 0), radius=central_radius, color=color.orange, emissive=True
    )
    equatorial_plane = ring(
        pos=vector(0, 0, 0),
        axis=vector(0, 0, 1),
        radius=apoapsis * 1.15,
        thickness=apoapsis * 0.002,
        color=color.gray(0.35),
    )
    polar_axis = cylinder(
        pos=vector(0, 0, -apoapsis * 0.3),
        axis=vector(0, 0, apoapsis * 0.6),
        radius=apoapsis * 0.002,
        color=color.gray(0.5),
    )
    return central_body, equatorial_plane, polar_axis


# --- Main driver --------------------------------------------------------

def main():
    params = gather_parameters()
    GM = params["GM"]
    r_periapsis = params["r_periapsis"]
    eccentricity = params["eccentricity"]
    precession_rad = math.radians(params["precession_deg"])
    inclination = math.radians(params["inclination_deg"])
    nodal_rate_rad = math.radians(params["nodal_rate_deg"])
    num_orbits = params["num_orbits"]

    a = semi_major_axis(r_periapsis, eccentricity)
    apoapsis = a * (1 + eccentricity)
    v0 = periapsis_speed(GM, r_periapsis, eccentricity)
    L = r_periapsis * v0
    c_squared = c_squared_for_target_precession(GM, a, eccentricity, precession_rad)
    period = orbital_period(GM, a)
    dt = period / N_STEPS_PER_ORBIT

    total_steps = N_STEPS_PER_ORBIT * num_orbits
    desired_frames = max(1, int(TARGET_ANIMATION_SECONDS * RENDER_FPS))
    steps_per_frame = max(1, round(total_steps / desired_frames))

    scene = setup_scene(apoapsis)
    central_radius = max(apoapsis * 0.02, 4.0)
    draw_reference_geometry(apoapsis, central_radius)

    raan0 = 0.0
    scene.append_to_caption(
        "Gray ring: reference plane.  Faint curve: unperturbed (non-precessing) "
        "Kepler ellipse.  Bright curve: actual simulated (precessing) trajectory.\n"
    )
    curve(pos=reference_ellipse_curve(a, eccentricity, inclination, raan0), color=color.gray(0.6), radius=central_radius * 0.1)

    satellite_radius = max(apoapsis * 0.012, 2.5)
    satellite = sphere(pos=vector(r_periapsis, 0, 0), radius=satellite_radius, color=color.cyan)
    trail = curve(radius=central_radius * 0.12, color=color.cyan)
    trail_n_points = {"value": 0}

    def trail_append(pos):
        trail.append(pos=pos)
        trail_n_points["value"] += 1
        if trail_n_points["value"] > MAX_TRAIL_POINTS:
            trail.pop(0)
            trail_n_points["value"] -= 1

    def trail_clear():
        while trail_n_points["value"] > 0:
            trail.pop(0)
            trail_n_points["value"] -= 1

    node_indicator = cylinder(
        pos=vector(0, 0, 0),
        axis=vector(apoapsis * 1.1, 0, 0),
        radius=apoapsis * 0.003,
        color=color.yellow,
    )

    state = {"running": True, "trail_on": True}

    def toggle_pause(b):
        state["running"] = not state["running"]
        b.text = "Resume" if not state["running"] else "Pause"

    def toggle_trail(b):
        state["trail_on"] = not state["trail_on"]
        if not state["trail_on"]:
            trail_clear()
        b.text = "Show trail" if not state["trail_on"] else "Hide trail"

    button(text="Pause", bind=toggle_pause)
    scene.append_to_caption("  ")
    button(text="Hide trail", bind=toggle_trail)
    scene.append_to_caption("\n")
    status = wtext(text="")

    r, r_dot, theta = r_periapsis, 0.0, 0.0
    r_history = [r]
    theta_history = [theta]
    periapsis_thetas = []
    target_deg = params["precession_deg"]

    def measured_precession_deg():
        if len(periapsis_thetas) < 2:
            return None
        diffs = [
            math.degrees(periapsis_thetas[i + 1] - periapsis_thetas[i]) - 360.0
            for i in range(len(periapsis_thetas) - 1)
        ]
        return sum(diffs) / len(diffs)

    def parabola_vertex(t0, y0, t1, y1, t2, y2):
        denom = (t0 - t1) * (t0 - t2) * (t1 - t2)
        aa = (t2 * (y1 - y0) + t1 * (y0 - y2) + t0 * (y2 - y1)) / denom
        bb = (t2 * t2 * (y0 - y1) + t1 * t1 * (y2 - y0) + t0 * t0 * (y1 - y2)) / denom
        return -bb / (2 * aa)

    steps_done = 0
    while steps_done < total_steps:
        rate(RENDER_FPS)
        if not state["running"]:
            continue

        for _ in range(steps_per_frame):
            if steps_done >= total_steps:
                break
            r, r_dot, theta = step_state(r, r_dot, theta, dt, GM, L, c_squared)
            r_history.append(r)
            theta_history.append(theta)
            if len(r_history) >= 3 and r_history[-2] < r_history[-3] and r_history[-2] < r_history[-1]:
                t_star = parabola_vertex(
                    theta_history[-3], r_history[-3],
                    theta_history[-2], r_history[-2],
                    theta_history[-1], r_history[-1],
                )
                periapsis_thetas.append(t_star)
            steps_done += 1

        elapsed_orbits = theta / (2 * math.pi)
        raan = raan0 + nodal_rate_rad * elapsed_orbits

        x_pf, y_pf = r * math.cos(theta), r * math.sin(theta)
        pos = perifocal_to_inertial(x_pf, y_pf, inclination, raan)
        satellite.pos = pos
        if state["trail_on"]:
            trail_append(pos)

        node_dir = vector(math.cos(raan), math.sin(raan), 0)
        node_indicator.axis = node_dir * apoapsis * 1.1

        measured = measured_precession_deg()
        measured_text = f"{measured:.2f}" if measured is not None else "collecting data..."
        status.text = (
            f"Orbit {elapsed_orbits:5.2f} / {num_orbits}   |   "
            f"RAAN (node) = {math.degrees(raan) % 360:6.1f} deg   |   "
            f"apsidal precession: target {target_deg:.2f} deg/orbit, measured {measured_text} deg/orbit"
        )

    final_measured = measured_precession_deg()
    print("\n=== Simulation complete ===")
    print(f"Target apsidal precession per orbit:   {target_deg:.4f} deg")
    if final_measured is not None:
        print(f"Measured apsidal precession per orbit: {final_measured:.4f} deg")
    else:
        print("Not enough periapsis passages were recorded to measure precession.")
    print("The 3D window remains open and interactive (drag to rotate, buttons still work).")

    while True:
        rate(30)


if __name__ == "__main__":
    main()
