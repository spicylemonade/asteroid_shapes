"""
Target Selection Query Engine

Identifies candidate asteroids for shape modeling based on:
1. NEO status (q < 1.3 AU) or large MBA (D > 100 km)
2. Sufficient lightcurve data in ALCDEF
3. Not already modeled in DAMIT
4. Priority scoring
"""

import csv
import gzip
import math
import os
import re
import zipfile
from collections import defaultdict

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Asteroids with models in the DAMIT database
DAMIT_MODELED = set()


def load_damit_list(damit_csv_path=None):
    """Load list of asteroids with models in DAMIT database.

    Parses the asteroids.csv table from a DAMIT database export.
    Falls back to a minimal hardcoded list if no export is available.
    """
    global DAMIT_MODELED

    # Try to load from DAMIT export CSV
    if damit_csv_path is None:
        damit_csv_path = os.path.join(REPO_ROOT, 'data', 'damit_asteroids.csv')

    if os.path.exists(damit_csv_path):
        DAMIT_MODELED = set()
        with open(damit_csv_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header:
                # Find the 'number' column index
                try:
                    num_idx = header.index('number')
                except ValueError:
                    num_idx = 1  # default: second column
                for row in reader:
                    try:
                        num = int(row[num_idx])
                        DAMIT_MODELED.add(num)
                    except (ValueError, IndexError):
                        continue
        return DAMIT_MODELED

    # Fallback: warn and return empty set
    print("WARNING: No DAMIT asteroids.csv found. DAMIT crossmatch disabled.")
    print(f"  Expected path: {damit_csv_path}")
    print("  Extract from DAMIT export: tar -xzf damit-*.tar.gz */tables/asteroids.csv")
    DAMIT_MODELED = set()
    return DAMIT_MODELED


def estimate_diameter(H, albedo=0.15):
    """Estimate diameter from absolute magnitude using D = 1329 / sqrt(pv) * 10^(-H/5)."""
    return 1329.0 / math.sqrt(albedo) * 10**(-H / 5.0)


def parse_alcdef_index(zip_path):
    """Quick index of ALCDEF data: count files and datapoints per asteroid number."""
    index = defaultdict(lambda: {'count': 0, 'datapoints': 0, 'name': ''})

    with zipfile.ZipFile(zip_path, 'r') as zf:
        for fname in zf.namelist():
            if not fname.endswith('.txt'):
                continue
            # Extract asteroid number from filename
            match = re.match(r'ALCDEF_(\d+)_(.+)\.txt', fname)
            if not match:
                continue
            num = int(match.group(1))
            name = match.group(2).replace('_', ' ')

            # Quick datapoint count
            try:
                content = zf.read(fname).decode('utf-8', errors='replace')
                dp_count = content.count('DATA=')
                sess_count = content.count('STARTMETADATA')
            except Exception:
                dp_count = 0
                sess_count = 0

            index[num]['count'] += sess_count
            index[num]['datapoints'] += dp_count
            if not index[num]['name']:
                index[num]['name'] = name

    return dict(index)


def select_candidates(mpcorb_path=None, alcdef_path=None,
                      min_lightcurves=2, min_datapoints=20,
                      max_candidates=100, verbose=True):
    """Select candidate asteroids for shape modeling.

    Returns list of candidate dicts sorted by priority score.
    """
    if mpcorb_path is None:
        mpcorb_path = os.path.join(REPO_ROOT, 'MPCORB.DAT.gz')
    if alcdef_path is None:
        alcdef_path = os.path.join(REPO_ROOT, 'ALCDEF_ALL.zip')

    load_damit_list()

    if verbose:
        print("Indexing ALCDEF data...")
    alcdef_index = parse_alcdef_index(alcdef_path)
    if verbose:
        print(f"  {len(alcdef_index)} asteroids with ALCDEF data")

    if verbose:
        print("Parsing MPCORB...")
    candidates = []

    with gzip.open(mpcorb_path, 'rt', errors='replace') as f:
        header_done = False
        for line in f:
            if line.strip().startswith('---'):
                header_done = True
                continue
            if not header_done or len(line) < 103:
                continue

            try:
                desig = line[0:7].strip()
                num = int(desig)
            except ValueError:
                continue

            try:
                H = float(line[8:13].strip())
                e = float(line[70:79].strip())
                a = float(line[92:103].strip())
            except (ValueError, IndexError):
                continue

            q = a * (1.0 - e)  # perihelion distance
            is_neo = q < 1.3
            diameter = estimate_diameter(H)
            is_large = diameter > 100.0

            # Priority 1: NEO or large MBA
            if not (is_neo or is_large):
                continue

            # Priority 3: Not in DAMIT
            if num in DAMIT_MODELED:
                continue

            # Priority 4: Has sufficient ALCDEF data
            if num not in alcdef_index:
                continue

            alc = alcdef_index[num]
            if alc['datapoints'] < min_datapoints:
                continue

            # Extract name
            name = alc['name'] if alc['name'] else ''
            if not name and len(line) > 166:
                nm = line[166:194].strip()
                m = re.search(r'\((\d+)\)\s*(.*)', nm)
                name = m.group(2).strip() if m and m.group(2).strip() else nm

            # Priority score: higher is better
            score = 0.0
            if is_neo:
                score += 10.0
            if is_large:
                score += 5.0
            score += min(alc['datapoints'] / 100.0, 5.0)  # data richness
            score += min(alc['count'] / 5.0, 3.0)  # session diversity

            candidates.append({
                'asteroid_id': num,
                'name': name,
                'neo_flag': is_neo,
                'estimated_diameter_km': round(diameter, 1),
                'H': H,
                'q_AU': round(q, 4),
                'a_AU': round(a, 4),
                'num_lightcurves': alc['count'],
                'num_datapoints': alc['datapoints'],
                'priority_score': round(score, 2),
            })

    # Sort by priority score descending
    candidates.sort(key=lambda x: x['priority_score'], reverse=True)
    candidates = candidates[:max_candidates]

    if verbose:
        print(f"Selected {len(candidates)} candidates")
        print(f"  NEOs: {sum(1 for c in candidates if c['neo_flag'])}")
        print(f"  Large MBAs: {sum(1 for c in candidates if not c['neo_flag'])}")

    return candidates


def save_candidates_csv(candidates, output_path):
    """Save candidates to CSV file."""
    fieldnames = ['asteroid_id', 'name', 'neo_flag', 'estimated_diameter_km',
                  'num_lightcurves', 'num_datapoints', 'priority_score']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(candidates)


if __name__ == '__main__':
    candidates = select_candidates(verbose=True)
    out_path = os.path.join(REPO_ROOT, 'results', 'candidate_list.csv')
    save_candidates_csv(candidates, out_path)
    print(f"\nSaved {len(candidates)} candidates to {out_path}")
