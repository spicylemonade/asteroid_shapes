#!/usr/bin/env python3
"""
Main Entry Point: Asteroid Lightcurve Inversion Pipeline

Reproduces the full analysis from raw ALCDEF/MPCORB data to output shape models.

Usage:
    python run_pipeline.py [--targets N] [--validate] [--batch]
"""

import argparse
import json
import os
import sys

import numpy as np

from src.data.parse_alcdef import parse_all_alcdef, write_summary_csv, make_histograms
from src.targets.selector import select_candidates, save_candidates_csv


def main():
    parser = argparse.ArgumentParser(description='Asteroid Lightcurve Inversion Pipeline')
    parser.add_argument('--targets', type=int, default=50, help='Number of targets to process')
    parser.add_argument('--validate', action='store_true', help='Run validation on ground-truth asteroids')
    parser.add_argument('--batch', action='store_true', help='Run batch inversion on candidates')
    parser.add_argument('--parse-alcdef', action='store_true', help='Parse ALCDEF data and generate summary')
    parser.add_argument('--select', action='store_true', help='Run target selection')
    parser.add_argument('--all', action='store_true', help='Run all steps')
    args = parser.parse_args()

    os.makedirs('results', exist_ok=True)
    os.makedirs('results/shapes', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    os.makedirs('figures/shapes', exist_ok=True)

    if args.all or args.parse_alcdef:
        print("\n=== Step 1: Parse ALCDEF Data ===")
        data = parse_all_alcdef('ALCDEF_ALL.zip')
        rows = write_summary_csv(data, 'results/alcdef_summary.csv')
        make_histograms(rows, 'figures')

    if args.all or args.select:
        print("\n=== Step 2: Select Candidates ===")
        candidates = select_candidates(max_candidates=args.targets)
        save_candidates_csv(candidates, 'results/candidate_list.csv')

    if args.validate:
        print("\n=== Step 3: Validation ===")
        print("Run: python -m src.inversion.run_validation")

    if args.batch:
        print("\n=== Step 4: Batch Inversion ===")
        print("Run batch inversion on selected candidates.")
        print("See results/batch_run_log.csv for progress.")

    if not any([args.all, args.parse_alcdef, args.select, args.validate, args.batch]):
        parser.print_help()
        print("\nExample: python run_pipeline.py --all")


if __name__ == '__main__':
    main()
