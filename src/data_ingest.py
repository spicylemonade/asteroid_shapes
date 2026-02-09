"""
ALCDEF data ingestion and lightcurve preprocessing module.
Reads ALCDEF text files, parses lightcurve data, computes viewing geometry
using MPCORB orbital elements.

References:
  - Kaasalainen et al. (2001) [Kaasalainen2001a in sources.bib]
  - Durech et al. (2009) [Durech2009 in sources.bib]
"""
import zipfile
import gzip
import os
import math
import numpy as np
from collections import defaultdict

# ============================================================
# ALCDEF Parsing
# ============================================================

def parse_alcdef_file(content):
    """Parse a single ALCDEF text file into lightcurve blocks.

    Returns list of dicts, each with 'meta' (dict) and 'data' (list of [jd, mag, err]).
    """
    blocks = []
    current_meta = {}
    data_points = []
    in_metadata = False

    for line in content.split('\n'):
        line = line.strip()
        if line == 'STARTMETADATA':
            if current_meta and data_points:
                blocks.append({'meta': current_meta, 'data': data_points})
            in_metadata = True
            current_meta = {}
            data_points = []
        elif line == 'ENDMETADATA':
            in_metadata = False
        elif in_metadata and '=' in line:
            key, _, val = line.partition('=')
            current_meta[key.strip()] = val.strip()
        elif line.startswith('DATA='):
            parts = line[5:].split('|')
            if len(parts) >= 2:
                try:
                    jd = float(parts[0])
                    mag = float(parts[1])
                    err = float(parts[2]) if len(parts) >= 3 else 0.05
                    data_points.append([jd, mag, err])
                except ValueError:
                    pass

    if current_meta and data_points:
        blocks.append({'meta': current_meta, 'data': data_points})

    return blocks


def load_alcdef_for_asteroid(zip_path, asteroid_number):
    """Load all lightcurve blocks for a specific numbered asteroid from ALCDEF ZIP.

    Parameters
    ----------
    zip_path : str
        Path to ALCDEF_ALL.zip
    asteroid_number : int or str
        Asteroid number (e.g. 433 for Eros)

    Returns
    -------
    list of dict
        Each dict has 'meta' and 'data' (Nx3 array: JD, mag, uncertainty)
    """
    asteroid_number = str(asteroid_number)
    zf = zipfile.ZipFile(zip_path)

    # Find matching file(s)
    prefix = f'ALCDEF_{asteroid_number}_'
    matching = [n for n in zf.namelist() if n.startswith(prefix)]

    if not matching:
        raise ValueError(f"No ALCDEF file found for asteroid {asteroid_number}")

    all_blocks = []
    for fname in matching:
        content = zf.read(fname).decode('utf-8', errors='replace')
        blocks = parse_alcdef_file(content)
        # Convert data to numpy arrays
        for block in blocks:
            block['data'] = np.array(block['data'], dtype=np.float64)
        all_blocks.extend(blocks)

    return all_blocks


# ============================================================
# MPCORB Orbital Elements
# ============================================================

def load_mpcorb_for_asteroid(gz_path, asteroid_number):
    """Load orbital elements for a specific asteroid from MPCORB.DAT.gz.

    Parameters
    ----------
    gz_path : str
        Path to MPCORB.DAT.gz
    asteroid_number : int or str
        Asteroid number

    Returns
    -------
    dict with orbital elements: a, e, i, node, peri, M, epoch, H, G
    """
    target = str(asteroid_number).zfill(5)
    header_passed = False

    with gzip.open(gz_path, 'rt', errors='replace') as f:
        for line in f:
            if not header_passed:
                if line.startswith('------'):
                    header_passed = True
                continue

            if len(line) < 103:
                continue

            designation = line[0:7].strip()
            # Match by number
            if designation == target:
                return _parse_mpcorb_record(line)

    raise ValueError(f"Asteroid {asteroid_number} not found in MPCORB")


def _parse_mpcorb_record(line):
    """Parse a single MPCORB line into orbital elements dict."""
    return {
        'designation': line[0:7].strip(),
        'H': _safe_float(line[8:13]),
        'G': _safe_float(line[14:19]),
        'epoch': line[20:25].strip(),
        'M_deg': _safe_float(line[26:35]),       # Mean anomaly
        'peri_deg': _safe_float(line[37:46]),     # Argument of perihelion
        'node_deg': _safe_float(line[48:57]),     # Longitude of ascending node
        'incl_deg': _safe_float(line[59:68]),     # Inclination
        'e': _safe_float(line[70:79]),            # Eccentricity
        'n_deg_day': _safe_float(line[80:91]),    # Mean daily motion
        'a_au': _safe_float(line[92:103]),        # Semimajor axis
        'name': line[166:194].strip() if len(line) > 194 else '',
    }


def _safe_float(s):
    try:
        return float(s.strip())
    except (ValueError, AttributeError):
        return None


# ============================================================
# Viewing Geometry Computation
# ============================================================

def kepler_solve(M_rad, e, tol=1e-10, max_iter=100):
    """Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly E."""
    E = M_rad
    for _ in range(max_iter):
        dE = (M_rad - E + e * math.sin(E)) / (1 - e * math.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E


def orbital_position_ecliptic(orb, jd):
    """Compute heliocentric ecliptic coordinates of an asteroid at given JD.

    Parameters
    ----------
    orb : dict
        Orbital elements from MPCORB (a_au, e, incl_deg, node_deg, peri_deg, M_deg, n_deg_day, epoch)
    jd : float
        Julian Date

    Returns
    -------
    np.array of shape (3,) - heliocentric ecliptic [x, y, z] in AU
    """
    # Unpack packed MPC epoch to approximate JD
    epoch_jd = unpack_epoch_to_jd(orb['epoch'])

    a = orb['a_au']
    e = orb['e']
    i = math.radians(orb['incl_deg'])
    node = math.radians(orb['node_deg'])
    peri = math.radians(orb['peri_deg'])
    n = orb['n_deg_day']  # mean motion deg/day

    # Mean anomaly at observation time
    dt = jd - epoch_jd
    M = math.radians(orb['M_deg'] + n * dt) % (2 * math.pi)

    # Solve Kepler's equation
    E = kepler_solve(M, e)

    # True anomaly
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2),
                         math.sqrt(1 - e) * math.cos(E / 2))

    # Distance
    r = a * (1 - e * math.cos(E))

    # Position in orbital plane
    x_orb = r * math.cos(nu)
    y_orb = r * math.sin(nu)

    # Rotation to ecliptic frame
    cos_node = math.cos(node)
    sin_node = math.sin(node)
    cos_peri = math.cos(peri)
    sin_peri = math.sin(peri)
    cos_i = math.cos(i)
    sin_i = math.sin(i)

    x = (cos_node * cos_peri - sin_node * sin_peri * cos_i) * x_orb + \
        (-cos_node * sin_peri - sin_node * cos_peri * cos_i) * y_orb
    y = (sin_node * cos_peri + cos_node * sin_peri * cos_i) * x_orb + \
        (-sin_node * sin_peri + cos_node * cos_peri * cos_i) * y_orb
    z = (sin_peri * sin_i) * x_orb + (cos_peri * sin_i) * y_orb

    return np.array([x, y, z])


def earth_position_ecliptic(jd):
    """Approximate Earth's heliocentric ecliptic position at given JD.
    Uses simplified orbital elements for Earth (sufficient for geometry computation).
    """
    # Earth orbital elements (J2000, approximate)
    a = 1.00000261  # AU
    e = 0.01671123
    i = math.radians(0.00005)
    node = math.radians(-11.26064)
    peri = math.radians(102.94719)
    L0 = 100.46435  # mean longitude at J2000 epoch
    n = 0.9856474  # mean motion deg/day

    # J2000 epoch
    epoch_jd = 2451545.0
    dt = jd - epoch_jd

    M = math.radians((L0 + n * dt - 102.94719) % 360.0)
    E = kepler_solve(M, e)

    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2),
                         math.sqrt(1 - e) * math.cos(E / 2))
    r = a * (1 - e * math.cos(E))

    x_orb = r * math.cos(nu)
    y_orb = r * math.sin(nu)

    cos_node = math.cos(node)
    sin_node = math.sin(node)
    cos_peri = math.cos(peri)
    sin_peri = math.sin(peri)
    cos_i = math.cos(i)
    sin_i = math.sin(i)

    x = (cos_node * cos_peri - sin_node * sin_peri * cos_i) * x_orb + \
        (-cos_node * sin_peri - sin_node * cos_peri * cos_i) * y_orb
    y = (sin_node * cos_peri + cos_node * sin_peri * cos_i) * x_orb + \
        (-sin_node * sin_peri + cos_node * cos_peri * cos_i) * y_orb
    z = (sin_peri * sin_i) * x_orb + (cos_peri * sin_i) * y_orb

    return np.array([x, y, z])


def compute_viewing_geometry(asteroid_pos, earth_pos):
    """Compute phase angle, aspect angle, and solar elongation.

    Parameters
    ----------
    asteroid_pos : np.array (3,) heliocentric ecliptic
    earth_pos : np.array (3,) heliocentric ecliptic

    Returns
    -------
    dict with phase_angle_deg, solar_elongation_deg, helio_dist_au, geo_dist_au
    """
    # Vectors
    sun_to_ast = asteroid_pos
    earth_to_ast = asteroid_pos - earth_pos
    sun_to_earth = earth_pos

    r_helio = np.linalg.norm(sun_to_ast)
    r_geo = np.linalg.norm(earth_to_ast)
    r_earth = np.linalg.norm(sun_to_earth)

    # Phase angle: angle Sun-Asteroid-Earth
    if r_helio > 0 and r_geo > 0:
        cos_phase = np.dot(-sun_to_ast, earth_to_ast) / (r_helio * r_geo)
        cos_phase = np.clip(cos_phase, -1, 1)
        phase_angle = math.degrees(math.acos(cos_phase))
    else:
        phase_angle = 0.0

    # Solar elongation: angle Sun-Earth-Asteroid
    if r_earth > 0 and r_geo > 0:
        cos_elong = np.dot(-sun_to_earth, earth_to_ast) / (r_earth * r_geo)
        cos_elong = np.clip(cos_elong, -1, 1)
        solar_elongation = math.degrees(math.acos(cos_elong))
    else:
        solar_elongation = 0.0

    return {
        'phase_angle_deg': phase_angle,
        'solar_elongation_deg': solar_elongation,
        'helio_dist_au': r_helio,
        'geo_dist_au': r_geo,
    }


def unpack_epoch_to_jd(packed):
    """Convert MPC packed epoch format to Julian Date.
    Format: e.g., K25BL
    Century: I=18, J=19, K=20
    Year: two digits
    Month: 1-9 or A=10, B=11, C=12
    Day: 1-9, A=10, B=11, ..., V=31
    """
    if not packed or len(packed) < 5:
        return 2451545.0  # Default J2000

    century_map = {'I': 18, 'J': 19, 'K': 20}
    month_map = {str(i): i for i in range(1, 10)}
    month_map.update({chr(65 + i - 10): i for i in range(10, 13)})  # A=10, B=11, C=12
    day_map = {str(i): i for i in range(1, 10)}
    day_map.update({chr(65 + i - 10): i for i in range(10, 32)})  # A=10 ... V=31

    try:
        century = century_map.get(packed[0], 20)
        year = century * 100 + int(packed[1:3])
        month = month_map.get(packed[3], 1)
        day = day_map.get(packed[4], 1)

        # Convert to JD (simplified)
        a = (14 - month) // 12
        y = year + 4800 - a
        m = month + 12 * a - 3
        jd = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
        return float(jd)
    except (KeyError, ValueError, IndexError):
        return 2451545.0


# ============================================================
# High-level preprocessing
# ============================================================

def preprocess_asteroid(zip_path, gz_path, asteroid_number):
    """Load and preprocess lightcurve data for a single asteroid.

    Returns
    -------
    dict with keys:
      'asteroid_number': int
      'blocks': list of preprocessed lightcurve block dicts
      'orbital_elements': dict
    Each block dict has:
      'meta': original ALCDEF metadata
      'data': Nx3 array [JD, mag, err]
      'geometry': list of dicts with viewing geometry per data point
      'reduced_mag': array of reduced magnitudes
    """
    # Load lightcurve data
    blocks = load_alcdef_for_asteroid(zip_path, asteroid_number)
    if not blocks:
        raise ValueError(f"No lightcurve data for asteroid {asteroid_number}")

    # Load orbital elements
    orb = load_mpcorb_for_asteroid(gz_path, asteroid_number)

    # Compute geometry and reduced magnitudes for each block
    for block in blocks:
        geom_list = []
        reduced_mags = []

        for row in block['data']:
            jd = row[0]
            mag = row[1]

            # Compute asteroid and Earth positions
            ast_pos = orbital_position_ecliptic(orb, jd)
            earth_pos = earth_position_ecliptic(jd)

            geom = compute_viewing_geometry(ast_pos, earth_pos)
            geom_list.append(geom)

            # Reduced magnitude: correct for distance
            r = geom['helio_dist_au']
            delta = geom['geo_dist_au']
            if r > 0 and delta > 0:
                reduced_mag = mag - 5 * math.log10(r * delta)
            else:
                reduced_mag = mag

            reduced_mags.append(reduced_mag)

        block['geometry'] = geom_list
        block['reduced_mag'] = np.array(reduced_mags)

    return {
        'asteroid_number': int(asteroid_number),
        'name': orb.get('name', ''),
        'blocks': blocks,
        'orbital_elements': orb,
    }


# ============================================================
# Test
# ============================================================

if __name__ == '__main__':
    import sys
    repo_root = os.path.dirname(os.path.dirname(__file__))
    zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
    gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')

    test_asteroids = [433, 1580, 1036]  # Eros, Betulia, Ganymed

    for ast_num in test_asteroids:
        print(f"\n{'='*60}")
        print(f"Processing asteroid {ast_num}")
        print(f"{'='*60}")
        try:
            result = preprocess_asteroid(zip_path, gz_path, ast_num)
            print(f"Name: {result['name']}")
            print(f"Orbital elements: a={result['orbital_elements']['a_au']:.4f} AU, "
                  f"e={result['orbital_elements']['e']:.4f}, "
                  f"i={result['orbital_elements']['incl_deg']:.2f} deg")
            print(f"Number of lightcurve blocks: {len(result['blocks'])}")
            total_pts = sum(len(b['data']) for b in result['blocks'])
            print(f"Total data points: {total_pts}")

            # Show geometry for first block
            if result['blocks']:
                b = result['blocks'][0]
                g = b['geometry'][0]
                print(f"First block session: {b['meta'].get('SESSIONDATE', '?')}")
                print(f"First point geometry: phase={g['phase_angle_deg']:.2f} deg, "
                      f"r={g['helio_dist_au']:.3f} AU, "
                      f"delta={g['geo_dist_au']:.3f} AU, "
                      f"elong={g['solar_elongation_deg']:.2f} deg")
                print(f"Phase angle from metadata: {b['meta'].get('PHASE', '?')} deg")
        except Exception as e:
            print(f"Error: {e}")
