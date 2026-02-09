"""
Parse ALCDEF_ALL.zip to build a comprehensive catalog of asteroid lightcurve data.
Output: results/alcdef_catalog.csv
"""
import zipfile
import csv
import os
import sys
from collections import defaultdict

def parse_alcdef_file(content):
    """Parse a single ALCDEF text file into lightcurve blocks."""
    blocks = []
    current_meta = {}
    data_points = []
    in_metadata = False

    for line in content.split('\n'):
        line = line.strip()
        if line == 'STARTMETADATA':
            # Save previous block if it has data
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
                data_points.append(parts)

    # Last block
    if current_meta and data_points:
        blocks.append({'meta': current_meta, 'data': data_points})

    return blocks


def build_catalog(zip_path, output_csv):
    """Parse all ALCDEF files and build catalog."""
    zf = zipfile.ZipFile(zip_path)
    names = sorted(zf.namelist())

    # Aggregate per asteroid
    asteroids = defaultdict(lambda: {
        'name': '',
        'num_lightcurves': 0,
        'total_data_points': 0,
        'jd_min': float('inf'),
        'jd_max': float('-inf'),
        'observers': set(),
        'filters': set(),
        'phases': [],
        'mpc_codes': set(),
    })

    total_blocks = 0
    processed = 0

    for fname in names:
        if not fname.endswith('.txt'):
            continue
        try:
            content = zf.read(fname).decode('utf-8', errors='replace')
        except Exception:
            continue

        blocks = parse_alcdef_file(content)

        for block in blocks:
            meta = block['meta']
            data = block['data']

            obj_num = meta.get('OBJECTNUMBER', '0')
            obj_name = meta.get('OBJECTNAME', '')

            key = (obj_num, obj_name)
            rec = asteroids[key]
            rec['name'] = obj_name
            rec['num_lightcurves'] += 1
            rec['total_data_points'] += len(data)
            total_blocks += 1

            # Parse JD range
            for dp in data:
                try:
                    jd = float(dp[0])
                    if jd < rec['jd_min']:
                        rec['jd_min'] = jd
                    if jd > rec['jd_max']:
                        rec['jd_max'] = jd
                except (ValueError, IndexError):
                    pass

            # Observer
            obs = meta.get('OBSERVERS', meta.get('CONTACTNAME', ''))
            if obs:
                rec['observers'].add(obs)

            # Filter
            filt = meta.get('FILTER', '')
            if filt:
                rec['filters'].add(filt)

            # MPC code
            mpc = meta.get('MPCCODE', '')
            if mpc:
                rec['mpc_codes'].add(mpc)

            # Phase angle
            phase = meta.get('PHASE', '')
            if phase:
                try:
                    rec['phases'].append(float(phase))
                except ValueError:
                    pass

        processed += 1
        if processed % 5000 == 0:
            print(f"  Processed {processed}/{len(names)} files, {total_blocks} blocks so far...")

    print(f"Total files processed: {processed}")
    print(f"Total lightcurve blocks: {total_blocks}")
    print(f"Unique asteroid keys: {len(asteroids)}")

    # Write CSV
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'asteroid_number', 'asteroid_name', 'number_of_lightcurves',
            'total_data_points', 'jd_min', 'jd_max', 'date_range',
            'observer_codes', 'filter_bands', 'num_apparitions'
        ])

        for (obj_num, obj_name), rec in sorted(asteroids.items(), key=lambda x: (-x[1]['num_lightcurves'], x[0])):
            jd_min = rec['jd_min'] if rec['jd_min'] != float('inf') else ''
            jd_max = rec['jd_max'] if rec['jd_max'] != float('-inf') else ''
            date_range = f"{jd_min}-{jd_max}" if jd_min and jd_max else ''

            # Estimate apparitions from JD spread (rough: group by ~365 day gaps)
            phases_sorted = sorted(rec['phases']) if rec['phases'] else []
            n_apparitions = estimate_apparitions(rec)

            writer.writerow([
                obj_num,
                obj_name,
                rec['num_lightcurves'],
                rec['total_data_points'],
                jd_min,
                jd_max,
                date_range,
                ';'.join(sorted(rec['observers'])),
                ';'.join(sorted(rec['filters'])),
                n_apparitions
            ])

    # Summary stats
    many_lc = [(k, v) for k, v in asteroids.items() if v['num_lightcurves'] >= 20]
    print(f"\nAsteroids with >=20 lightcurves: {len(many_lc)}")
    for (num, name), v in sorted(many_lc, key=lambda x: -x[1]['num_lightcurves'])[:20]:
        print(f"  ({num}) {name}: {v['num_lightcurves']} lightcurves, {v['total_data_points']} points")

    return asteroids


def estimate_apparitions(rec):
    """Estimate number of distinct apparitions from JD timestamps.
    Group observations separated by >120 days as separate apparitions."""
    if rec['jd_min'] == float('inf'):
        return 0
    jd_range = rec['jd_max'] - rec['jd_min']
    if jd_range < 120:
        return 1
    # Rough estimate
    return max(1, int(jd_range / 365.25) + 1)


if __name__ == '__main__':
    zip_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'ALCDEF_ALL.zip')
    output_csv = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results', 'alcdef_catalog.csv')
    print(f"Parsing {zip_path}...")
    build_catalog(zip_path, output_csv)
    print(f"\nCatalog saved to {output_csv}")
