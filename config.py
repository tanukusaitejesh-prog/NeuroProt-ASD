"""
Central Configuration for ASD Gene Prioritization Pipeline
"""

import os
from pathlib import Path
import torch

BASE_DIR = Path(__file__).resolve().parent

# Output directories
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
CACHE_DIR = BASE_DIR / "cache"

for d in [RESULTS_DIR, FIGURES_DIR, CACHE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Hardware device
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Random seed for reproducibility
RANDOM_SEED = 42

# ESM-2 Models
# Lightweight model for fast local verification / development
ESM2_TINY_MODEL = "facebook/esm2_t6_8M_UR50D"
# Full model for production runs
ESM2_FULL_MODEL = "facebook/esm2_t33_650M_UR50D"

# Target disorders for multi-task modeling (§5.3)
TARGET_DISORDERS = ["ASD", "ID", "SCZ", "EPI"]

# Conformal prediction parameters (§6)
CONFORMAL_TARGET_COVERAGE = 0.90
CONFORMAL_CALIBRATION_FRACTION = 0.15

# Top-N threshold for forecASD comparison (§7.3)
TOP_N_CANDIDATES = 500
