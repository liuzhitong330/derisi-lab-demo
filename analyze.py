#!/usr/bin/env python3
"""Descriptive PhIP-seq control-reference audit; no clinical classification.

Requires pandas/numpy, two CC0 Dryad human-proteome CSVs, and a locally fetched
phagepy checkout at the pinned commit. Only original wrappers are distributed;
phagepy has no explicit redistribution license in the inspected checkout.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import statistics
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

COMMIT = "d600a664c69541f8b04fe47fdce166fe52c16f6f"
GENES = ["SNX8", "ERFL", "KDELR1", "UBE3A"]
ALIASES = {"LOC390937": "ERFL"}  # NCBI Gene 390937; do not confuse with ERF.
HASHES = {
    "Human_Proteome_MIS-C_patients.csv": "599682c2eba758b91d88a8d15c147feba7553637025922730108ec15f6e53f4c",
    "Human_Proteome_at-risk_controls.csv": "bc9de7b03040c0017f28d5890f6a4c22a1e5fb8588d1fd70c571bc7d0e8baab8",
}

class MatrixAdapter:
    """Minimal documented AnnData-like adapter for this pure numeric function.

    Does not run preprocessing, read raw counts, or pretend FC values are RPK.
    The unchanged upstream function only uses obs, var, layers, obsm, and row
    subsetting; these are provided explicitly here to avoid unrelated dependencies.
    """
    def __init__(self, matrix, ids, groups, genes):
        self.obs = pd.DataFrame({"group": groups}, index=ids)
        self.var = pd.DataFrame(index=genes)
        self.layers = {"FC_over_mock_IP": np.asarray(matrix, dtype=float)}
        self.obsm = {}

    def __getitem__(self, mask):
        indices = np.asarray(mask, dtype=bool)
        return MatrixAdapter(self.layers["FC_over_mock_IP"][indices],
            self.obs.index[indices], self.obs.loc[indices, "group"], self.var.index)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, default=Path(__file__).parent / "inputs")
    parser.add_argument("--phagepy", type=Path, default=Path(__file__).parent / "phagepy")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "derived")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    actual_commit = subprocess.check_output(["git", "-C", str(args.phagepy), "rev-parse", "HEAD"], text=True).strip()
    assert actual_commit == COMMIT, (actual_commit, COMMIT)
    assert hashlib.sha256((args.phagepy / "phagepy" / "analysis.py").read_bytes()).hexdigest() == "29400ca115cc8a78743f209a1f8df92bd2da501b0c49322af04df84c0aaf6396"
    for filename, expected in HASHES.items():
        assert hashlib.sha256((args.inputs / filename).read_bytes()).hexdigest() == expected, filename
    frames = [pd.read_csv(args.inputs / filename, index_col=0) for filename in HASHES]
    cases, controls = frames
    audit = {}
    for label, frame in zip(["MIS-C", "at-risk control"], frames):
        assert frame.index.is_unique and frame.columns.is_unique
        assert np.isfinite(frame.to_numpy()).all() and (frame.to_numpy() >= 0).all()
        audit[label] = {"genes": len(frame), "samples": len(frame.columns), "missing": int(frame.isna().sum().sum())}
    assert cases.index.equals(controls.index), "Gene order mismatch"
    assert not set(cases.columns) & set(controls.columns), "Group ID overlap"
    assert cases.shape[1] == 199 and controls.shape[1] == 45
    assert "ERFL" not in cases.index and "LOC390937" in cases.index
    cases = cases.rename(index=ALIASES)
    controls = controls.rename(index=ALIASES)
    assert set(GENES).issubset(cases.index)
    # Execute unchanged, inspected upstream analysis module. No raw-data normalization.
    spec = importlib.util.spec_from_file_location("phagepy_analysis", args.phagepy / "phagepy" / "analysis.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    matrix = pd.concat([cases.loc[GENES].T, controls.loc[GENES].T])
    adapter = MatrixAdapter(matrix.to_numpy(), matrix.index,
        ["MIS-C"] * 199 + ["at-risk control"] * 45, GENES)
    module.generate_SS_hits(adapter, "group", "at-risk control", layer="FC_over_mock_IP", z_cutoff=3)
    ctl = controls.loc[GENES].T.to_numpy()
    assert (ctl.std(axis=0, ddof=0) > 0).all()
    expected_z = (matrix.to_numpy() - ctl.mean(axis=0)) / ctl.std(axis=0, ddof=0)
    assert np.allclose(expected_z, adapter.layers["Z_at-risk control"], rtol=1e-12, atol=1e-12)
    assert np.array_equal(expected_z > 3, adapter.obsm["SS_at-risk control_Z-3"].to_numpy())
    records, summaries = [], []
    data = {"meta": {"dataset": "doi:10.7272/Q6SJ1HVH", "dataset_version": 293273,
        "repository": "https://github.com/h-s-miller/phagepy", "commit": COMMIT,
        "input_unit": "gene-level fold change over mock-IP, after peptide-read aggregation and read-percentage normalization", "control_sd_ddof": 0,
        "cutoffs": [2, 3, 4, 5], "primary_genes": GENES[:3], "comparison_gene": "UBE3A",
        "cases": 199, "controls": 45, "data_audit": audit,
        "gene_aliases": {"LOC390937": {"display": "ERFL", "source": "https://www.ncbi.nlm.nih.gov/gene/390937"}},
        "interpretation": "Above-reference is a descriptive screen, not antibody confirmation, disease prediction, or a clinical positive.",
        "missing_covariates": ["IVIG status", "collection timing", "technical replicates", "raw counts", "linked orthogonal assay", "material availability"]}, "genes": []}
    for gene in GENES:
        x = cases.loc[gene].to_numpy()
        y = controls.loc[gene].to_numpy()
        mu, sd = float(y.mean()), float(y.std(ddof=0))
        assert abs(mu - statistics.mean(y.tolist())) < 1e-12
        assert abs(sd - statistics.pstdev(y.tolist())) < 1e-12
        loo = np.array([[np.delete(y, j).mean(), np.delete(y, j).std(ddof=0)] for j in range(len(y))])
        gene_data = {"gene": gene, "mean": mu, "sd": sd,
            "cases": [{"id": sid, "fc": float(v)} for sid, v in zip(cases.columns, x)],
            "controls": [{"id": sid, "fc": float(v), "heldout_z": float((v - loo[j, 0]) / loo[j, 1])} for j, (sid, v) in enumerate(zip(controls.columns, y))],
            "leave_one_out": [{"omitted_id": sid, "mean": float(a), "sd": float(b)} for sid, (a, b) in zip(controls.columns, loo)],
            "settings": {}}
        for k in [2, 3, 4, 5]:
            threshold = mu + k * sd
            ts = loo[:, 0] + k * loo[:, 1]
            hit = x[:, None] > ts[None, :]
            hits = hit.sum(axis=1)
            stable_above = hits == 45
            sensitive = (hits > 0) & (hits < 45)
            stable_below = hits == 0
            heldout_control_hit = y > ts
            row = {"gene": gene, "k": k, "threshold": threshold, "loo_min": float(ts.min()), "loo_max": float(ts.max()),
                "full_reference_above": int((x > threshold).sum()), "stable_above": int(stable_above.sum()),
                "sensitive": int(sensitive.sum()), "stable_below": int(stable_below.sum()),
                "heldout_controls_above": int(heldout_control_hit.sum()),
                "largest_cutoff_shift_id": str(controls.columns[np.argmax(np.abs(ts - threshold))]),
                "largest_cutoff_shift": float(ts[np.argmax(np.abs(ts - threshold))] - threshold)}
            summaries.append(row)
            gene_data["settings"][str(k)] = row
            for j, (sid, value) in enumerate(zip(cases.columns, x)):
                records.append({"gene": gene, "k": k, "sample_id": sid, "cohort": "MIS-C", "fold_change_over_mock_IP": value,
                    "full_reference_z": (value - mu) / sd, "loo_above_count": int(hits[j]), "loo_runs": 45,
                    "reference_category": "stable_above" if stable_above[j] else "sensitive" if sensitive[j] else "stable_below",
                    "status": "Proposed assay planning only; IVIG, collection timing and material availability unverified"})
        data["genes"].append(gene_data)
    # Minimal exact numeric independent check of upstream result, one value per gene.
    checks = [{"gene": gene, "sample": cases.columns[0], "standard_library_z":
        (float(cases.loc[gene].iloc[0]) - statistics.mean(controls.loc[gene].tolist())) / statistics.pstdev(controls.loc[gene].tolist()),
        "phagepy_z": float(adapter.layers["Z_at-risk control"][0, j])} for j, gene in enumerate(GENES)]
    assert all(abs(c["standard_library_z"] - c["phagepy_z"]) < 1e-12 for c in checks)
    data["verification"] = {"unchanged_upstream_function_executed": "generate_SS_hits",
        "adapter": "Minimal AnnData-compatible numeric adapter, not the upstream raw-count loading pipeline",
        "full_matrix_z_parity": True, "full_matrix_hit_parity": True, "spot_checks": checks}
    (args.out / "demo-data.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    pd.DataFrame(records).to_csv(args.out / "sample_reference_audit.csv", index=False)
    pd.DataFrame(summaries).to_csv(args.out / "threshold_summary.csv", index=False)
    (args.out / "verification.json").write_text(json.dumps(data["verification"], indent=2, allow_nan=False) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False))
    print(json.dumps(audit))

if __name__ == "__main__":
    main()
