"""
Parse MPCORB.DAT.gz to extract orbital elements for all minor planets.
Flags NEOs and estimates diameters from H magnitude.
Output: results/mpcorb_parsed.csv
"""
import gzip
import csv
import os
import math


def unpack_mpc_epoch(packed):
    """Convert packed MPC epoch format to JD-like string.
    Format: K25BL = 2025-11-21 (century, year, month, day in packed form)
    We just store the packed format as-is for now."""
    return packed.strip()


def estimate_diameter_km(H, albedo=0.15):
    """Estimate asteroid diameter from absolute magnitude H and assumed albedo.
    D(km) = 1329 / sqrt(albedo) * 10^(-H/5)
    """
    if H is None or H == '':
        return None
    try:
        H = float(H)
        return 1329.0 / math.sqrt(albedo) * (10 ** (-H / 5.0))
    except (ValueError, ZeroDivisionError):
        return None


def is_neo(a, e, q=None):
    """Classify as NEO based on orbital elements.
    NEO: perihelion q < 1.3 AU (q = a*(1-e))
    """
    try:
        a = float(a)
        e = float(e)
        q_val = a * (1 - e) if q is None else float(q)
        return q_val < 1.3
    except (ValueError, TypeError):
        return False


def classify_orbit(a, e):
    """Classify orbit type: NEO, MBA, or other."""
    try:
        a = float(a)
        e = float(e)
        q = a * (1 - e)
        Q = a * (1 + e)

        if q < 1.3:
            return 'NEO'
        elif 2.0 <= a <= 3.5 and e < 0.4:
            return 'MBA'
        elif a > 30:
            return 'TNO'
        else:
            return 'OTHER'
    except (ValueError, TypeError):
        return 'UNKNOWN'


def parse_mpcorb_line(line):
    """Parse a single MPCORB.DAT fixed-width line.
    See: https://www.minorplanetcenter.org/iau/info/MPOrbitFormat.html
    """
    if len(line) < 103:
        return None

    try:
        designation = line[0:7].strip()
        H_str = line[8:13].strip()
        G_str = line[14:19].strip()
        epoch = line[20:25].strip()
        M = line[26:35].strip()       # Mean anomaly
        peri = line[37:46].strip()     # Argument of perihelion
        node = line[48:57].strip()     # Longitude of ascending node
        incl = line[59:68].strip()     # Inclination
        e_str = line[70:79].strip()    # Eccentricity
        n = line[80:91].strip()        # Mean daily motion
        a_str = line[92:103].strip()   # Semimajor axis

        # Extended fields (may be shorter)
        name = line[166:194].strip() if len(line) > 194 else ''

        H = float(H_str) if H_str else None
        G = float(G_str) if G_str else None
        a = float(a_str) if a_str else None
        e = float(e_str) if e_str else None

        neo_flag = is_neo(a, e) if a and e else False
        orbit_type = classify_orbit(a, e) if a and e else 'UNKNOWN'
        est_diameter = estimate_diameter_km(H) if H else None

        return {
            'designation': designation,
            'name': name,
            'H': H,
            'G': G,
            'epoch': epoch,
            'mean_anomaly': M,
            'arg_perihelion': peri,
            'long_asc_node': node,
            'inclination': incl,
            'eccentricity': e_str,
            'mean_daily_motion': n,
            'semimajor_axis': a_str,
            'neo_flag': neo_flag,
            'orbit_type': orbit_type,
            'est_diameter_km': round(est_diameter, 2) if est_diameter else None,
            'large_flag': est_diameter is not None and est_diameter > 100,
        }
    except Exception:
        return None


def parse_mpcorb(gz_path, output_csv):
    """Parse entire MPCORB.DAT.gz file."""
    records = []
    header_passed = False
    count = 0
    neo_count = 0
    large_count = 0

    with gzip.open(gz_path, 'rt', errors='replace') as f:
        for line in f:
            if not header_passed:
                if line.startswith('------'):
                    header_passed = True
                continue

            if len(line.strip()) == 0:
                continue

            rec = parse_mpcorb_line(line)
            if rec:
                records.append(rec)
                if rec['neo_flag']:
                    neo_count += 1
                if rec['large_flag']:
                    large_count += 1
                count += 1

            if count % 200000 == 0:
                print(f"  Parsed {count} records...")

    print(f"Total records parsed: {count}")
    print(f"NEOs: {neo_count}")
    print(f"Diameter > 100km: {large_count}")

    # Write CSV
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    fieldnames = [
        'designation', 'name', 'H', 'G', 'epoch', 'mean_anomaly',
        'arg_perihelion', 'long_asc_node', 'inclination', 'eccentricity',
        'mean_daily_motion', 'semimajor_axis', 'neo_flag', 'orbit_type',
        'est_diameter_km', 'large_flag'
    ]

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Saved to {output_csv}")

    # Show some NEO examples
    neos = [r for r in records if r['neo_flag']][:10]
    print(f"\nSample NEOs:")
    for r in neos:
        print(f"  {r['designation']} {r['name']}: H={r['H']}, a={r['semimajor_axis']}, D~{r['est_diameter_km']}km")

    # Show large asteroids
    larges = sorted([r for r in records if r['large_flag']], key=lambda x: -(x['est_diameter_km'] or 0))[:15]
    print(f"\nLargest asteroids (D>100km):")
    for r in larges:
        print(f"  {r['designation']} {r['name']}: H={r['H']}, D~{r['est_diameter_km']}km")

    return records


if __name__ == '__main__':
    gz_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'MPCORB.DAT.gz')
    output_csv = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results', 'mpcorb_parsed.csv')
    parse_mpcorb(gz_path, output_csv)
