"""
MPCORB Parser and Ephemeris Geometry Calculator

Parses MPCORB.DAT.gz orbital elements and computes viewing geometry
for asteroid observations using two-body Keplerian propagation.
"""

import gzip
import math
import os
import re
from datetime import datetime

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Constants
GM_SUN = 1.32712440018e20  # m^3/s^2
AU = 1.495978707e11  # meters
DAY = 86400.0  # seconds
GM_SUN_AU = 0.01720209895**2  # AU^3/day^2 (Gaussian gravitational constant squared)
DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


def unpack_mpc_epoch(packed):
    """Unpack MPC packed epoch format to JD.

    Format: CYYHM where C=century letter, YY=year, H=half-month, M=day
    Century: I=18, J=19, K=20
    Half-month: 1-9,A-C (Jan-Dec)
    Day: 1-9,A-V (1-31)
    """
    if len(packed) != 5:
        return 2451545.0  # J2000 default

    century_map = {'I': 1800, 'J': 1900, 'K': 2000}
    month_map = {str(i): i for i in range(1, 10)}
    month_map.update({'A': 10, 'B': 11, 'C': 12})
    day_map = {str(i): i for i in range(1, 10)}
    for i, c in enumerate('ABCDEFGHIJKLMNOPQRSTUV'):
        day_map[c] = 10 + i

    try:
        year = century_map.get(packed[0], 2000) + int(packed[1:3])
        month = month_map.get(packed[3], 1)
        day = day_map.get(packed[4], 1)
        return calendar_to_jd(year, month, day)
    except (KeyError, ValueError):
        return 2451545.0


def calendar_to_jd(year, month, day):
    """Convert calendar date to Julian Date."""
    if month <= 2:
        year -= 1
        month += 12
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    return int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5


def parse_mpcorb(filepath=None):
    """Parse MPCORB.DAT.gz into a dictionary keyed by asteroid number and designation."""
    if filepath is None:
        filepath = os.path.join(REPO_ROOT, 'MPCORB.DAT.gz')

    asteroids = {}
    header_done = False

    open_func = gzip.open if filepath.endswith('.gz') else open

    with open_func(filepath, 'rt', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('---'):
                header_done = True
                continue
            if not header_done or len(line) < 103:
                continue

            try:
                desig = line[0:7].strip()
                H = float(line[8:13].strip()) if line[8:13].strip() else 99.0
                G = float(line[14:19].strip()) if line[14:19].strip() else 0.15
                epoch_packed = line[20:25].strip()
                M = float(line[26:35].strip()) if line[26:35].strip() else 0.0
                peri = float(line[37:46].strip()) if line[37:46].strip() else 0.0
                node = float(line[48:57].strip()) if line[48:57].strip() else 0.0
                incl = float(line[59:68].strip()) if line[59:68].strip() else 0.0
                e = float(line[70:79].strip()) if line[70:79].strip() else 0.0
                n = float(line[80:91].strip()) if line[80:91].strip() else 0.0
                a = float(line[92:103].strip()) if line[92:103].strip() else 0.0
            except (ValueError, IndexError):
                continue

            # Extract name from columns 166+
            name = ''
            if len(line) > 166:
                name = line[166:194].strip()
                # Clean up parenthetical name
                m = re.search(r'\((\d+)\)\s*(.*)', name)
                if m:
                    name = m.group(2).strip() if m.group(2).strip() else name

            epoch_jd = unpack_mpc_epoch(epoch_packed)
            # Perihelion distance
            q = a * (1.0 - e)

            record = {
                'designation': desig,
                'H': H,
                'G': G,
                'epoch_jd': epoch_jd,
                'M0': M * DEG2RAD,
                'peri': peri * DEG2RAD,
                'node': node * DEG2RAD,
                'incl': incl * DEG2RAD,
                'e': e,
                'n': n * DEG2RAD,  # mean motion in rad/day
                'a': a,
                'q': q,
                'name': name,
            }

            # Key by number if it's a numbered asteroid
            try:
                num = int(desig)
                asteroids[num] = record
                asteroids[desig] = record
            except ValueError:
                asteroids[desig] = record

            if name:
                asteroids[name] = record

    return asteroids


def solve_kepler(M, e, tol=1e-12, max_iter=100):
    """Solve Kepler's equation M = E - e*sin(E) for E."""
    # Initial guess
    E = M + e * np.sin(M)
    for _ in range(max_iter):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E += dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def kepler_position(a, e, incl, node, peri, M):
    """Compute heliocentric ecliptic position from orbital elements.

    All angles in radians. Returns position in AU (ecliptic J2000).
    """
    E = solve_kepler(M, e)
    # True anomaly
    nu = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(E / 2.0),
        np.sqrt(1.0 - e) * np.cos(E / 2.0)
    )
    # Heliocentric distance
    r = a * (1.0 - e * np.cos(E))

    # Position in orbital plane
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)

    # Rotation to ecliptic frame
    cos_peri = np.cos(peri)
    sin_peri = np.sin(peri)
    cos_node = np.cos(node)
    sin_node = np.sin(node)
    cos_incl = np.cos(incl)
    sin_incl = np.sin(incl)

    x_ecl = (cos_node * cos_peri - sin_node * sin_peri * cos_incl) * x_orb + \
             (-cos_node * sin_peri - sin_node * cos_peri * cos_incl) * y_orb
    y_ecl = (sin_node * cos_peri + cos_node * sin_peri * cos_incl) * x_orb + \
             (-sin_node * sin_peri + cos_node * cos_peri * cos_incl) * y_orb
    z_ecl = (sin_peri * sin_incl) * x_orb + (cos_peri * sin_incl) * y_orb

    return np.array([x_ecl, y_ecl, z_ecl])


def earth_position(jd):
    """Approximate Earth heliocentric ecliptic position at given JD.

    Uses simplified Keplerian orbit for Earth.
    """
    # Earth orbital elements (J2000 mean, approximate)
    a_earth = 1.00000261  # AU
    e_earth = 0.01671123
    incl_earth = 0.00005 * DEG2RAD  # nearly zero
    node_earth = -11.26064 * DEG2RAD
    peri_earth = 102.93768 * DEG2RAD
    L0_earth = 100.46457 * DEG2RAD  # mean longitude at J2000

    # J2000 epoch
    jd_j2000 = 2451545.0
    dt = jd - jd_j2000  # days

    # Mean motion
    n_earth = 2.0 * np.pi / 365.25636  # rad/day

    # Mean longitude
    L = L0_earth + n_earth * dt
    # Mean anomaly
    M = L - peri_earth
    M = M % (2.0 * np.pi)

    return kepler_position(a_earth, e_earth, incl_earth, node_earth,
                           peri_earth - node_earth, M)


def compute_geometry(asteroid_record, obs_jd):
    """Compute observation geometry for an asteroid at given epoch.

    Parameters
    ----------
    asteroid_record : dict
        Orbital elements from parse_mpcorb()
    obs_jd : float
        Julian Date of observation

    Returns
    -------
    dict with keys:
        r_ast : heliocentric position of asteroid (AU, ecliptic)
        r_earth : heliocentric position of Earth (AU, ecliptic)
        delta : geocentric vector to asteroid (AU)
        r_helio : heliocentric distance (AU)
        delta_geo : geocentric distance (AU)
        phase_angle : solar phase angle (degrees)
        solar_elongation : solar elongation angle (degrees)
        sun_dir_ecl : unit vector toward Sun in ecliptic frame
        obs_dir_ecl : unit vector toward observer in ecliptic frame
    """
    rec = asteroid_record
    dt = obs_jd - rec['epoch_jd']
    M = rec['M0'] + rec['n'] * dt
    M = M % (2.0 * np.pi)

    r_ast = kepler_position(rec['a'], rec['e'], rec['incl'],
                            rec['node'], rec['peri'], M)
    r_earth = earth_position(obs_jd)

    delta = r_ast - r_earth  # geocentric vector
    r_helio = np.linalg.norm(r_ast)
    delta_geo = np.linalg.norm(delta)

    # Phase angle: angle Sun-Asteroid-Observer
    # cos(phase) = (-r_ast . delta) / (|r_ast| * |delta|)
    cos_phase = np.dot(-r_ast, delta) / (r_helio * delta_geo)
    cos_phase = np.clip(cos_phase, -1.0, 1.0)
    phase_angle = np.arccos(cos_phase) * RAD2DEG

    # Solar elongation: angle Sun-Earth-Asteroid
    cos_elong = np.dot(-r_earth, delta) / (np.linalg.norm(r_earth) * delta_geo)
    cos_elong = np.clip(cos_elong, -1.0, 1.0)
    solar_elongation = np.arccos(cos_elong) * RAD2DEG

    # Direction vectors (unit) in ecliptic frame
    sun_dir = -r_ast / r_helio  # from asteroid toward Sun
    obs_dir = -delta / delta_geo  # from asteroid toward observer (Earth)

    return {
        'r_ast': r_ast,
        'r_earth': r_earth,
        'delta': delta,
        'r_helio': r_helio,
        'delta_geo': delta_geo,
        'phase_angle': phase_angle,
        'solar_elongation': solar_elongation,
        'sun_dir_ecl': sun_dir,
        'obs_dir_ecl': obs_dir,
    }


def compute_aspect_angles(sun_dir_ecl, obs_dir_ecl, pole_lambda, pole_beta):
    """Compute sub-observer and sub-solar latitude/longitude on asteroid.

    Parameters
    ----------
    sun_dir_ecl : array, unit vector toward Sun in ecliptic
    obs_dir_ecl : array, unit vector toward observer in ecliptic
    pole_lambda : float, pole ecliptic longitude in degrees
    pole_beta : float, pole ecliptic latitude in degrees

    Returns
    -------
    dict with sub_obs_lat, sub_obs_lon, sub_sun_lat, sub_sun_lon (degrees)
    """
    lam = pole_lambda * DEG2RAD
    bet = pole_beta * DEG2RAD

    # Spin axis in ecliptic
    sz = np.array([np.cos(bet) * np.cos(lam),
                   np.cos(bet) * np.sin(lam),
                   np.sin(bet)])

    # Sub-observer latitude = angle between observer dir and equatorial plane
    sub_obs_lat = np.arcsin(np.clip(np.dot(obs_dir_ecl, sz), -1, 1)) * RAD2DEG

    # For longitude, we need a reference direction in the equatorial plane
    # Use x-axis projection
    sx = np.array([-np.sin(lam), np.cos(lam), 0.0])
    sy = np.cross(sz, sx)
    sy = sy / np.linalg.norm(sy)
    sx = np.cross(sy, sz)

    obs_x = np.dot(obs_dir_ecl, sx)
    obs_y = np.dot(obs_dir_ecl, sy)
    sub_obs_lon = np.arctan2(obs_y, obs_x) * RAD2DEG

    sub_sun_lat = np.arcsin(np.clip(np.dot(sun_dir_ecl, sz), -1, 1)) * RAD2DEG
    sun_x = np.dot(sun_dir_ecl, sx)
    sun_y = np.dot(sun_dir_ecl, sy)
    sub_sun_lon = np.arctan2(sun_y, sun_x) * RAD2DEG

    return {
        'sub_obs_lat': sub_obs_lat,
        'sub_obs_lon': sub_obs_lon % 360,
        'sub_sun_lat': sub_sun_lat,
        'sub_sun_lon': sub_sun_lon % 360,
    }
