"""
ALCDEF Data Parser

Parses all ALCDEF lightcurve files from ALCDEF_ALL.zip and produces
summary statistics and histogram plots.
"""

import csv
import os
import re
import sys
import zipfile
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

try:
    import seaborn as sns
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
except ImportError:
    pass

mpl.rcParams.update({
    'figure.figsize': (8, 5),
    'figure.dpi': 300,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
    'font.family': 'serif',
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def parse_alcdef_file(content):
    """Parse a single ALCDEF file content into list of session records."""
    sessions = []
    current_session = None
    in_metadata = False

    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue

        if line == 'STARTMETADATA':
            in_metadata = True
            current_session = {
                'object_number': None,
                'object_name': '',
                'session_date': '',
                'filter_band': '',
                'observer': '',
                'datapoints': [],
            }
            continue

        if line == 'ENDMETADATA':
            in_metadata = False
            continue

        if in_metadata and '=' in line:
            key, _, value = line.partition('=')
            key = key.strip()
            value = value.strip()

            if key == 'OBJECTNUMBER':
                try:
                    current_session['object_number'] = int(value)
                except ValueError:
                    current_session['object_number'] = 0
            elif key == 'OBJECTNAME':
                current_session['object_name'] = value
            elif key == 'SESSIONDATE':
                current_session['session_date'] = value
            elif key == 'FILTER':
                current_session['filter_band'] = value
            elif key == 'OBSERVERS':
                current_session['observer'] = value
            elif key == 'MAGBAND':
                if not current_session['filter_band']:
                    current_session['filter_band'] = value

        elif line.startswith('DATA='):
            parts = line[5:].split('|')
            if len(parts) >= 2:
                try:
                    jd = float(parts[0])
                    mag = float(parts[1])
                    err = float(parts[2]) if len(parts) > 2 else 0.0
                    if current_session is not None:
                        current_session['datapoints'].append((jd, mag, err))
                except (ValueError, IndexError):
                    pass

        # When we hit next STARTMETADATA or EOF, save current session
        if line == 'STARTMETADATA' and current_session is not None and current_session.get('datapoints'):
            sessions.append(current_session)
            # Don't reset - it was already reset above

    # Save last session
    if current_session is not None and current_session.get('datapoints'):
        sessions.append(current_session)

    return sessions


def parse_all_alcdef(zip_path):
    """Parse all ALCDEF files from zip archive."""
    # Keyed by (object_number, object_name)
    asteroid_data = defaultdict(lambda: {
        'sessions': [],
        'all_jds': [],
        'filters': set(),
        'observers': set(),
        'total_datapoints': 0,
    })

    with zipfile.ZipFile(zip_path, 'r') as zf:
        file_list = [n for n in zf.namelist() if n.endswith('.txt')]
        total = len(file_list)
        print(f"Parsing {total} ALCDEF files...")

        for idx, fname in enumerate(file_list):
            if (idx + 1) % 5000 == 0:
                print(f"  Progress: {idx + 1}/{total} ({100*(idx+1)/total:.1f}%)")

            try:
                content = zf.read(fname).decode('utf-8', errors='replace')
            except Exception as e:
                print(f"  Warning: Could not read {fname}: {e}")
                continue

            sessions = parse_alcdef_file(content)

            for sess in sessions:
                obj_num = sess['object_number'] if sess['object_number'] else 0
                obj_name = sess['object_name']
                key = (obj_num, obj_name)

                rec = asteroid_data[key]
                rec['sessions'].append(sess)
                jds = [dp[0] for dp in sess['datapoints']]
                rec['all_jds'].extend(jds)
                rec['total_datapoints'] += len(sess['datapoints'])
                if sess['filter_band']:
                    rec['filters'].add(sess['filter_band'])
                if sess['observer']:
                    rec['observers'].add(sess['observer'])

    print(f"Parsed {total} files. Found {len(asteroid_data)} unique asteroids.")
    return asteroid_data


def write_summary_csv(asteroid_data, output_path):
    """Write summary CSV file."""
    rows = []
    for (obj_num, obj_name), rec in sorted(asteroid_data.items()):
        jds = rec['all_jds']
        if jds:
            min_jd = min(jds)
            max_jd = max(jds)
            # Convert JD to approximate calendar date
            date_range = f"JD {min_jd:.1f} - {max_jd:.1f}"
        else:
            date_range = ""

        rows.append({
            'asteroid_number': obj_num,
            'name': obj_name,
            'num_sessions': len(rec['sessions']),
            'total_datapoints': rec['total_datapoints'],
            'date_range': date_range,
            'filter_bands': ';'.join(sorted(rec['filters'])),
            'observers': ';'.join(sorted(rec['observers'])),
        })

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'asteroid_number', 'name', 'num_sessions', 'total_datapoints',
            'date_range', 'filter_bands', 'observers'
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Summary CSV written to {output_path} ({len(rows)} asteroids)")
    return rows


def make_histograms(rows, figures_dir):
    """Create publication-quality histogram plots."""
    try:
        palette = sns.color_palette("deep")
    except Exception:
        palette = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3']

    # Histogram 1: Datapoints per asteroid
    datapoints = [r['total_datapoints'] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    bins = np.logspace(0, np.log10(max(datapoints) + 1), 50)
    ax.hist(datapoints, bins=bins, color=palette[0], edgecolor='white',
            linewidth=0.5, alpha=0.85)
    ax.set_xscale('log')
    ax.set_xlabel('Total Data Points per Asteroid')
    ax.set_ylabel('Number of Asteroids')
    ax.set_title('Distribution of Photometric Data Points per Asteroid')
    ax.axvline(100, color=palette[3], linestyle='--', linewidth=1.5,
               label='100 points threshold')
    ax.legend(frameon=True, shadow=True)

    median_dp = np.median(datapoints)
    ax.annotate(f'Median: {median_dp:.0f}', xy=(median_dp, ax.get_ylim()[1] * 0.8),
                fontsize=10, color=palette[1], fontweight='bold')

    fig.savefig(os.path.join(figures_dir, 'alcdef_datapoints_histogram.png'), dpi=300)
    fig.savefig(os.path.join(figures_dir, 'alcdef_datapoints_histogram.pdf'))
    plt.close(fig)
    print(f"Saved datapoints histogram")

    # Histogram 2: Number of sessions per asteroid
    sessions = [r['num_sessions'] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    max_sess = min(max(sessions), 200)
    bins = np.arange(0, max_sess + 2, max(1, max_sess // 50))
    ax.hist(sessions, bins=bins, color=palette[1], edgecolor='white',
            linewidth=0.5, alpha=0.85)
    ax.set_xlabel('Number of Observation Sessions per Asteroid')
    ax.set_ylabel('Number of Asteroids')
    ax.set_title('Distribution of Observation Sessions per Asteroid')

    median_s = np.median(sessions)
    ax.annotate(f'Median: {median_s:.0f}', xy=(median_s * 1.5, ax.get_ylim()[1] * 0.8),
                fontsize=10, color=palette[0], fontweight='bold')

    fig.savefig(os.path.join(figures_dir, 'alcdef_sessions_histogram.png'), dpi=300)
    fig.savefig(os.path.join(figures_dir, 'alcdef_sessions_histogram.pdf'))
    plt.close(fig)
    print(f"Saved sessions histogram")

    # Histogram 3: Temporal coverage (date range span in years)
    spans = []
    for r in rows:
        dr = r['date_range']
        if dr and ' - ' in dr:
            parts = dr.split(' - ')
            try:
                jd1 = float(parts[0].replace('JD ', ''))
                jd2 = float(parts[1].replace('JD ', ''))
                span_years = (jd2 - jd1) / 365.25
                if span_years >= 0:
                    spans.append(span_years)
            except ValueError:
                pass

    if spans:
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        ax.hist(spans, bins=50, color=palette[2], edgecolor='white',
                linewidth=0.5, alpha=0.85)
        ax.set_xlabel('Temporal Coverage (years)')
        ax.set_ylabel('Number of Asteroids')
        ax.set_title('Distribution of Observation Time Span')

        median_span = np.median(spans)
        ax.annotate(f'Median: {median_span:.1f} yr', xy=(median_span * 1.2, ax.get_ylim()[1] * 0.8),
                    fontsize=10, color=palette[4], fontweight='bold')

        fig.savefig(os.path.join(figures_dir, 'alcdef_temporal_coverage.png'), dpi=300)
        fig.savefig(os.path.join(figures_dir, 'alcdef_temporal_coverage.pdf'))
        plt.close(fig)
        print(f"Saved temporal coverage histogram")

    # Summary statistics
    print(f"\n=== ALCDEF Summary Statistics ===")
    print(f"Total unique asteroids: {len(rows)}")
    print(f"Total data points: {sum(datapoints)}")
    print(f"Median data points/asteroid: {np.median(datapoints):.0f}")
    print(f"Mean data points/asteroid: {np.mean(datapoints):.0f}")
    print(f"Max data points: {max(datapoints)} ({rows[np.argmax(datapoints)]['name']})")
    print(f"Median sessions/asteroid: {np.median(sessions):.0f}")
    if spans:
        print(f"Median temporal span: {np.median(spans):.1f} years")


def main():
    zip_path = os.path.join(REPO_ROOT, 'ALCDEF_ALL.zip')
    results_dir = os.path.join(REPO_ROOT, 'results')
    figures_dir = os.path.join(REPO_ROOT, 'figures')

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    if not os.path.exists(zip_path):
        print(f"ERROR: {zip_path} not found")
        sys.exit(1)

    asteroid_data = parse_all_alcdef(zip_path)
    rows = write_summary_csv(asteroid_data, os.path.join(results_dir, 'alcdef_summary.csv'))
    make_histograms(rows, figures_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
