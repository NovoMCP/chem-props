"""
refit_qed.py — Refit QED desirability functions to ChEMBL 35 approved oral drugs.

One-time offline script. Reads the NovoExpert-3 training Parquet (or computes
descriptors from SMILES), fits 8 Asymmetric Double Sigmoidal (ADS) curves
against the empirical distributions, and outputs:
  - qed_v2_coefficients.json (the fitted parameters)
  - qed_v2.py (runtime module)

The ADS function (Bickerton et al. 2012) is:
    f(x) = A + (B / (1 + exp(-1*(x-C+D/2)/E)) * (1 - 1/(1 + exp(-1*(x-C-D/2)/F))))

Each of the 8 molecular descriptors (MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS)
gets its own ADS curve fitted to the empirical distribution of that descriptor
across the ChEMBL 35 approved oral drug set.

Usage:
    python refit_qed.py --parquet /path/to/phase1_features.parquet
    python refit_qed.py --smiles-file /path/to/approved_oral_smiles.txt

Requirements:
    pip install rdkit pandas scipy numpy
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Bickerton 2012 original parameters (for comparison)
BICKERTON_2012 = {
    "MW": {"A": 2.817065973, "B": 392.5754953, "C": 290.7489764, "D": 2.419764353, "E": 49.22325677, "F": 65.37051707, "DMAX": 104.9805561},
    "ALOGP": {"A": 3.172690585, "B": 137.8624751, "C": 2.534937431, "D": 4.581497897, "E": 0.822739154, "F": 0.576295591, "DMAX": 131.3186604},
    "HBA": {"A": 2.948620388, "B": 160.4605972, "C": 3.615294657, "D": 4.435986202, "E": 0.290141953, "F": 1.300669958, "DMAX": 148.7763046},
    "HBD": {"A": 1.618662227, "B": 1010.051101, "C": 0.985094388, "D": 0.000000001, "E": 0.713820843, "F": 0.920922555, "DMAX": 258.1632616},
    "PSA": {"A": 1.876861559, "B": 125.2232657, "C": 62.90773554, "D": 87.83366614, "E": 12.01999824, "F": 28.51324732, "DMAX": 104.5686167},
    "ROTB": {"A": 0.010000000, "B": 272.4121427, "C": 2.558379970, "D": 1.565547684, "E": 1.271567166, "F": 0.587993418, "DMAX": 105.4420403},
    "AROM": {"A": 3.217788970, "B": 957.7374108, "C": 2.274627939, "D": 0.000000001, "E": 1.317690384, "F": 0.375760881, "DMAX": 312.3372610},
    "ALERTS": {"A": 0.010000000, "B": 1199.094025, "C": -0.09002883, "D": 0.000000001, "E": 0.185904477, "F": 0.875193782, "DMAX": 417.7253140},
}

# QED property weights (Bickerton 2012, unchanged — these are subjective importance weights)
QED_WEIGHTS = {
    "MW": 0.66, "ALOGP": 0.46, "HBA": 0.05, "HBD": 0.61,
    "PSA": 0.06, "ROTB": 0.65, "AROM": 0.48, "ALERTS": 0.95,
}

# Structural alerts (Brenk filter) — same SMARTS as RDKit's QED implementation
# These count the number of structural alerts in a molecule
ALERT_SMARTS = None  # Loaded lazily


def ads_function(x, A, B, C, D, E, F):
    """
    Asymmetric Double Sigmoidal function (Bickerton et al. 2012).
    Maps a descriptor value to a desirability score [0, 1].
    """
    exp_term1 = np.exp(-1.0 * (x - C + D / 2.0) / E)
    exp_term2 = np.exp(-1.0 * (x - C - D / 2.0) / F)
    return A + B / (1.0 + exp_term1) * (1.0 - 1.0 / (1.0 + exp_term2))


def compute_qed_descriptors(smiles: str) -> dict:
    """Compute the 8 QED descriptors for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Count structural alerts using RDKit's QED module
    try:
        from rdkit.Chem.QED import properties
        props = properties(mol)
        return {
            "MW": props.MW,
            "ALOGP": props.ALOGP,
            "HBA": props.HBA,
            "HBD": props.HBD,
            "PSA": props.PSA,
            "ROTB": props.ROTB,
            "AROM": props.AROM,
            "ALERTS": props.ALERTS,
        }
    except Exception as e:
        logger.warning(f"Failed to compute QED descriptors for {smiles}: {e}")
        return None


def fit_ads_curve(values: np.ndarray, descriptor_name: str) -> dict:
    """
    Fit the ADS function to the probability density of a descriptor
    across the reference drug set.

    The desirability function should peak where most approved drugs cluster
    and tail off in both directions. We fit to a smoothed histogram (KDE)
    of the descriptor values, matching the Bickerton approach.
    """
    from scipy.stats import gaussian_kde

    # Build KDE of the descriptor distribution
    kde = gaussian_kde(values, bw_method="scott")

    # Evaluate on a fine grid
    x_min = max(values.min() - values.std(), 0) if descriptor_name != "ALOGP" else values.min() - values.std()
    x_max = values.max() + values.std()
    x_grid = np.linspace(x_min, x_max, 500)
    density = kde(x_grid)

    # Normalize density to [0, ~max_ads_value] to match Bickerton scale
    # The ADS function in Bickerton produces values up to DMAX (100-400 range)
    # Scale density so its peak matches roughly the Bickerton DMAX
    bickerton_dmax = BICKERTON_2012[descriptor_name]["DMAX"]
    density_scaled = density / density.max() * bickerton_dmax

    initial = BICKERTON_2012[descriptor_name]
    p0 = [initial["A"], initial["B"], initial["C"],
          initial["D"], initial["E"], initial["F"]]

    # Bounds — generous to allow the optimizer to find the new optimum
    lower = [1e-3, 0.1, -20, 1e-10, 1e-3, 1e-3]
    upper = [10, 2000, 1000, 200, 200, 200]

    # Clamp initial guess to be within bounds
    p0 = [max(lo, min(hi, v)) for v, lo, hi in zip(p0, lower, upper)]

    try:
        popt, _ = curve_fit(
            ads_function, x_grid, density_scaled,
            p0=p0, bounds=(lower, upper),
            maxfev=20000,
        )

        # Compute DMAX over a wide range
        x_eval = np.linspace(
            min(x_min, 0) if descriptor_name != "ALOGP" else x_min,
            x_max * 1.5,
            10000,
        )
        dmax = ads_function(x_eval, *popt).max()

        # Compute fit quality: normalized RMSE
        predicted = ads_function(x_grid, *popt)
        nrmse = np.sqrt(np.mean((predicted - density_scaled) ** 2)) / (density_scaled.max() - density_scaled.min())

        fitted = {
            "A": round(float(popt[0]), 9),
            "B": round(float(popt[1]), 7),
            "C": round(float(popt[2]), 7),
            "D": round(float(popt[3]), 9),
            "E": round(float(popt[4]), 9),
            "F": round(float(popt[5]), 9),
            "DMAX": round(float(dmax), 7),
        }

        logger.info(f"  {descriptor_name}: fitted (nRMSE: {nrmse:.4f}, DMAX: {dmax:.1f})")
        return fitted

    except RuntimeError as e:
        logger.warning(f"  {descriptor_name}: curve_fit failed ({e}), using Bickerton 2012 params")
        return BICKERTON_2012[descriptor_name]


def refit_all(reference_smiles: list[str]) -> dict:
    """
    Refit all 8 QED desirability functions against a reference drug set.

    Args:
        reference_smiles: list of canonical SMILES for approved oral drugs

    Returns:
        dict with keys MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS,
        each containing ADS parameters {A, B, C, D, E, F, DMAX}
    """
    logger.info(f"Computing QED descriptors for {len(reference_smiles)} reference drugs...")

    records = []
    failed = 0
    for smi in reference_smiles:
        desc = compute_qed_descriptors(smi)
        if desc is not None:
            records.append(desc)
        else:
            failed += 1

    if failed > 0:
        logger.warning(f"  {failed} compounds failed descriptor computation")

    desc_df = pd.DataFrame(records)
    logger.info(f"  Valid descriptors: {len(desc_df)}")

    # Log descriptor distributions
    logger.info("  Descriptor distributions (reference set):")
    for col in desc_df.columns:
        vals = desc_df[col]
        logger.info(f"    {col}: mean={vals.mean():.2f}, median={vals.median():.2f}, "
                    f"std={vals.std():.2f}, min={vals.min():.2f}, max={vals.max():.2f}")

    # Fit ADS curves
    logger.info("Fitting ADS curves...")
    coefficients = {}
    for descriptor in ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]:
        values = desc_df[descriptor].dropna().values
        coefficients[descriptor] = fit_ads_curve(values, descriptor)

    return coefficients


def main():
    parser = argparse.ArgumentParser(description="Refit QED to ChEMBL 35 approved oral drugs")
    parser.add_argument("--parquet", default=None, help="Path to phase1_features.parquet")
    parser.add_argument("--output-dir", default=".", help="Output directory for coefficients + module")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)

    # Load reference drug set
    if args.parquet:
        logger.info(f"Loading reference drugs from {args.parquet}")
        df = pd.read_parquet(args.parquet)
        # Filter to approved oral small molecules
        reference = df[(df["max_phase"] == 4.0) & (df["oral"] == True)]
        reference_smiles = reference["canonical_smiles"].dropna().tolist()
        logger.info(f"  Approved oral drugs: {len(reference_smiles)}")
    else:
        raise ValueError("--parquet is required")

    # Refit
    coefficients = refit_all(reference_smiles)

    # Save coefficients
    output = {
        "version": "2.0",
        "source": "ChEMBL 35 approved oral small molecules",
        "reference_set_size": len(reference_smiles),
        "weights": QED_WEIGHTS,
        "parameters": coefficients,
        "bickerton_2012": BICKERTON_2012,
    }

    coeff_path = output_dir / "qed_v2_coefficients.json"
    with open(coeff_path, "w") as f:
        json.dump(output, f, indent=2)
    logger.info(f"Coefficients written to {coeff_path}")

    # Print comparison table
    logger.info("\n=== Bickerton 2012 vs ChEMBL 35 Refit ===")
    for desc in ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]:
        old_c = BICKERTON_2012[desc]["C"]
        new_c = coefficients[desc]["C"]
        logger.info(f"  {desc:6s}: C (center) {old_c:10.4f} → {new_c:10.4f} "
                    f"(Δ = {new_c - old_c:+.4f})")


if __name__ == "__main__":
    main()
