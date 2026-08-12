"""
qed_v2.py — QED Quantitative Estimate of Drug-likeness (ChEMBL 35 Refit)

Drop-in replacement for RDKit's QED.qed(mol). Uses desirability function
coefficients refit to 1,883 approved oral small molecules from ChEMBL 35,
replacing the Bickerton 2012 coefficients (fit to ~1,200 drugs from ChEMBL 9).

Key changes vs Bickerton 2012:
  - MW center shifted 291 → 343 (modern drugs are heavier)
  - ALOGP center shifted 2.5 → 3.1 (slightly more lipophilic)
  - PSA center shifted 63 → 49 (wider PSA acceptance)
  - All 8 ADS curves refit to ChEMBL 35 probability densities

Usage:
    from qed_v2 import qed_v2
    score = qed_v2(mol)  # returns float 0-1, same interface as QED.qed()

Coefficients source: refit_qed.py against ChEMBL 35 approved oral drugs.
Training labels: CC BY-SA 3.0 (ChEMBL attribution required in model metadata).
"""

import math
from rdkit.Chem.QED import properties as qed_properties

# Asymmetric Double Sigmoidal parameters, refit to ChEMBL 35 (2025)
# Each tuple: (A, B, C, D, E, F, DMAX)
# Source: refit_qed.py fitting ADS to KDE of 1,883 approved oral drugs
# ALERTS uses Bickerton 2012 params (structural alert definitions unchanged,
# and the discrete 0-6 count distribution doesn't fit ADS well)
_ADS_PARAMS_V2 = {
    "MW":     (0.884114097, 172.0047935, 343.4310508, 199.999999999, 59.02264291, 99.653312603, 107.2148574),
    "ALOGP":  (0.664084475, 219.209599, 3.0974692, 3.03429839, 1.490305987, 1.010007432, 132.3687256),
    "HBA":    (1.08721076, 206.8518499, 3.9536783, 4.590345172, 0.700874744, 2.164305312, 151.7160363),
    "HBD":    (0.295416937, 817.7190038, 0.4697693, 0.0, 0.559968122, 1.521183568, 246.2523991),
    "PSA":    (0.235189056, 378.4868114, 48.7765338, 0.0, 22.824926403, 47.056239983, 105.8743413),
    "ROTB":   (0.145587114, 295.905369, 3.2867822, 1.306419212, 1.561980294, 2.832847523, 103.5012958),
    "AROM":   (1.253254379, 819.889652, 2.2099896, 0.0, 1.133107445, 0.665790621, 219.7875439),
    "ALERTS": (0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140),
}

# Property weights (unchanged from Bickerton 2012 — these are subjective)
_WEIGHTS = (0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)
_WEIGHT_SUM = sum(_WEIGHTS)

# Descriptor names in the order returned by rdkit.Chem.QED.properties()
_DESCRIPTOR_ORDER = ("MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS")


def _ads(x: float, A: float, B: float, C: float, D: float, E: float, F: float) -> float:
    """Asymmetric Double Sigmoidal function."""
    try:
        exp1 = math.exp(-1.0 * (x - C + D / 2.0) / E)
    except OverflowError:
        exp1 = float("inf")
    try:
        exp2 = math.exp(-1.0 * (x - C - D / 2.0) / F)
    except OverflowError:
        exp2 = float("inf")
    return A + B / (1.0 + exp1) * (1.0 - 1.0 / (1.0 + exp2))


def qed_v2(mol) -> float:
    """
    Compute QED v2 score for an RDKit Mol object.

    Drop-in replacement for rdkit.Chem.QED.qed(mol).
    Returns float in [0, 1].
    """
    props = qed_properties(mol)

    log_sum = 0.0
    for i, (name, value) in enumerate(zip(_DESCRIPTOR_ORDER, props)):
        params = _ADS_PARAMS_V2[name]
        dmax = params[6]
        raw = _ads(value, *params[:6])
        desirability = raw / dmax
        desirability = max(desirability, 1e-10)  # avoid log(0)
        log_sum += _WEIGHTS[i] * math.log(desirability)

    return math.exp(log_sum / _WEIGHT_SUM)
