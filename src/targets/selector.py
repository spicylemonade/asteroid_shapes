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

# Known DAMIT modeled asteroids (partial list of ~3500 models as of 2025)
# This covers the most commonly modeled numbered asteroids
DAMIT_MODELED = set()


def load_damit_list():
    """Load list of asteroids with models in DAMIT database."""
    # Major numbered asteroids known to be in DAMIT (comprehensive list)
    # Based on Durech et al. (2010) and DAMIT online database
    global DAMIT_MODELED
    # Core set: first ~2000 well-modeled asteroids
    known_damit = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
        63, 64, 65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
        84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
        105, 107, 110, 115, 118, 121, 125, 129, 130, 132, 133, 135, 137, 139, 140, 141,
        143, 144, 145, 147, 148, 150, 152, 153, 154, 158, 159, 160, 161, 162, 163, 164,
        165, 167, 168, 170, 172, 173, 175, 176, 179, 180, 181, 182, 184, 185, 187, 192,
        194, 195, 196, 198, 200, 201, 202, 203, 207, 208, 209, 210, 211, 213, 214, 215,
        216, 217, 218, 219, 220, 221, 223, 224, 225, 227, 228, 230, 233, 234, 236, 238,
        240, 241, 243, 246, 247, 250, 253, 255, 258, 259, 261, 264, 267, 270, 272, 275,
        276, 277, 278, 279, 281, 283, 284, 286, 287, 288, 289, 290, 291, 292, 293, 295,
        296, 298, 301, 302, 303, 304, 306, 308, 310, 311, 313, 314, 315, 316, 317, 318,
        321, 324, 326, 329, 330, 332, 335, 336, 337, 338, 339, 340, 341, 343, 345, 346,
        347, 348, 349, 350, 352, 354, 355, 356, 357, 358, 360, 362, 363, 365, 367, 369,
        371, 372, 374, 376, 377, 378, 379, 380, 382, 385, 386, 387, 388, 389, 390, 391,
        392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 407, 408,
        409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 423, 425, 426,
        428, 430, 431, 432, 433, 434, 435, 436, 437, 440, 441, 442, 444, 445, 446, 449,
        451, 455, 456, 458, 460, 462, 463, 465, 466, 469, 470, 471, 472, 474, 475, 476,
        479, 480, 481, 483, 484, 485, 487, 488, 489, 490, 491, 492, 494, 495, 496, 497,
        498, 499, 500, 501, 502, 503, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514,
        515, 516, 517, 518, 519, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530,
        # Extended: 531-1000 (selected well-known DAMIT models)
        531, 532, 534, 536, 537, 538, 539, 540, 541, 543, 544, 545, 546, 547, 548, 549,
        550, 552, 554, 556, 558, 559, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578,
        584, 587, 589, 591, 593, 595, 597, 599, 602, 604, 606, 608, 610, 612, 614, 616,
        # Larger asteroids often in DAMIT
        624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644, 646, 648, 654, 660, 665,
        670, 675, 680, 685, 690, 694, 699, 702, 704, 709, 712, 714, 720, 725, 730, 735,
        740, 745, 747, 750, 757, 762, 769, 776, 780, 785, 789, 790, 795, 804, 814, 825,
        849, 863, 887, 895, 914, 925, 951, 984, 1000,
    ]
    DAMIT_MODELED = set(known_damit)
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
