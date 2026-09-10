"""
External Validation Orchestrator (§7.3 forecASD, §7.4 IMPC, §7.5 MAGMA)
"""

import subprocess
import sys
from pathlib import Path
import pandas as pd

from config import BASE_DIR, RESULTS_DIR


def run_forecasd_benchmarking(
    our_scores_csv: str,
    top_n: int = 500
) -> pd.DataFrame:
    """Executes §7.3 forecASD comparison script."""
    script = BASE_DIR / "01_forecasd_comparison.py"
    cmd = [
        sys.executable,
        str(script),
        "--our-scores", str(our_scores_csv),
        "--top-n", str(top_n)
    ]
    print(f"[INFO] Executing §7.3 forecASD comparison...")
    res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", check=True)
    print(res.stdout)

    novel_path = RESULTS_DIR / "novel_candidates.csv"
    assert novel_path.exists(), f"Expected {novel_path} to be generated."
    return pd.read_csv(novel_path)


def run_impc_enrichment(
    novel_candidates_csv: str
) -> dict:
    """Executes §7.4 IMPC mouse phenotype enrichment script."""
    script = BASE_DIR / "02_impc_enrichment.py"
    cmd = [
        sys.executable,
        str(script),
        "--novel-candidates", str(novel_candidates_csv)
    ]
    print(f"[INFO] Executing §7.4 IMPC mouse phenotype enrichment...")
    res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", check=True)
    print(res.stdout)

    summary_path = RESULTS_DIR / "impc_enrichment_summary.json"
    assert summary_path.exists(), f"Expected {summary_path} to be generated."
    import json
    with open(summary_path, "r") as f:
        return json.load(f)
