"""
Quick blind inversion test for 433 Eros using real ALCDEF data.

Loads ALCDEF photometry for Eros, runs convex lightcurve inversion starting
from the known period, compares recovered shape/pole against a ground-truth
mesh, and saves the result to results/.

Outputs:
    results/blind_test_433_Eros.obj  -- recovered shape mesh
    results/blind_test_433_Eros.json -- validation metrics and parameters
"""

import sys, json, time
sys.path.insert(0, '.')
import numpy as np
from lci_engine.parsers import load_alcdef_asteroid
from lci_engine.forward_model import create_triaxial_ellipsoid, load_mesh_obj, save_mesh_obj
from lci_engine.inversion import convex_inversion
from lci_engine.validation import compute_validation_report

# Load real ALCDEF data for Eros
sessions = load_alcdef_asteroid('ALCDEF_ALL.zip', asteroid_number=433)
print(f"Loaded {len(sessions)} sessions for 433 Eros")

# Known parameters for Eros
known = {
    'pole_lambda_deg': 11.35, 'pole_beta_deg': 17.22, 
    'period_hours': 5.27025,
    'a_km': 17.2, 'b_km': 5.6, 'c_km': 5.6
}

# Run blind convex inversion (use known period as starting point)
print("Running convex inversion...")
start = time.time()
result = convex_inversion(
    sessions,
    period_init=known['period_hours'],
    n_facets=200,
    lambda_smooth=0.3,
    max_iter=150,
    seed=42
)
elapsed = time.time() - start

print(f"Inversion completed in {elapsed:.1f}s")
print(f"  Period: {result.period:.5f} h (known: {known['period_hours']})")
print(f"  Pole: lambda={np.degrees(result.pole_lambda):.1f}°, beta={np.degrees(result.pole_beta):.1f}°")
print(f"  Residual RMS: {result.residual_rms:.4f} mag")

# Load ground truth
truth_verts, truth_faces = load_mesh_obj('results/ground_truth_433_Eros.obj')

# Compute validation metrics
report = compute_validation_report(
    result.vertices, result.faces, truth_verts, truth_faces,
    result.pole_lambda, result.pole_beta, result.period,
    np.radians(known['pole_lambda_deg']), np.radians(known['pole_beta_deg']),
    known['period_hours'], result.residual_rms
)

print(f"\nValidation Report:")
for key, val in report.items():
    print(f"  {key}: {val:.4f}")

# Save result
save_mesh_obj(result.vertices, result.faces, 'results/blind_test_433_Eros.obj')

validation_result = {
    'asteroid': '433 Eros',
    'recovered_period_h': result.period,
    'recovered_pole_lambda_deg': float(np.degrees(result.pole_lambda)),
    'recovered_pole_beta_deg': float(np.degrees(result.pole_beta)),
    'residual_rms_mag': float(result.residual_rms),
    'n_iterations': result.n_iterations,
    'converged': result.converged,
    'chi2': float(result.chi2),
    'elapsed_seconds': elapsed,
    'validation': report,
    'known': known,
}

with open('results/blind_test_433_Eros.json', 'w') as f:
    json.dump(validation_result, f, indent=2)

print("\nResults saved to results/blind_test_433_Eros.json")
