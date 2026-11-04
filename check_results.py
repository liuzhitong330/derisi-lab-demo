#!/usr/bin/env python3
"""Independent standard-library verification against the source CSVs."""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

def read_selected(path):
    values = {}
    with path.open(newline="") as handle:
        rows = csv.reader(handle)
        names = next(rows)[1:]
        assert len(names) == len(set(names))
        for row in rows:
            if row[0] in {"SNX8", "LOC390937", "KDELR1", "UBE3A"}:
                gene = "ERFL" if row[0] == "LOC390937" else row[0]
                values[gene] = dict(zip(names, map(float, row[1:])))
    return values

def close(a, b):
    assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-11), (a, b)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, default=Path(__file__).parent / "inputs")
    parser.add_argument("--derived", type=Path, default=Path(__file__).parent / "derived")
    args = parser.parse_args()
    data = json.loads((args.derived / "demo-data.json").read_text())
    cases = read_selected(args.inputs / "Human_Proteome_MIS-C_patients.csv")
    controls = read_selected(args.inputs / "Human_Proteome_at-risk_controls.csv")
    verified_scenarios = 0
    for marker in data["genes"]:
        gene = marker["gene"]
        case_values, control_values = cases[gene], controls[gene]
        assert len(case_values) == 199 and len(control_values) == 45
        for record in marker["cases"]:
            close(record["fc"], case_values[record["id"]])
        for record in marker["controls"]:
            close(record["fc"], control_values[record["id"]])
        ys = list(control_values.values())
        mu, sd = statistics.mean(ys), statistics.pstdev(ys)
        close(mu, marker["mean"])
        close(sd, marker["sd"])
        for k in (2, 3, 4, 5):
            row = marker["settings"][str(k)]
            cutoff = mu + k * sd
            omissions = {}
            for omitted_id in control_values:
                remainder = [v for sid, v in control_values.items() if sid != omitted_id]
                omissions[omitted_id] = statistics.mean(remainder) + k * statistics.pstdev(remainder)
            above = [sum(x > t for t in omissions.values()) for x in case_values.values()]
            expected = {
                "full_reference_above": sum(x > cutoff for x in case_values.values()),
                "stable_above": sum(n == 45 for n in above),
                "sensitive": sum(0 < n < 45 for n in above),
                "stable_below": sum(n == 0 for n in above),
                "heldout_controls_above": sum(v > omissions[sid] for sid, v in control_values.items()),
            }
            for name, value in expected.items():
                assert value == row[name], (gene, k, name, value, row[name])
            close(cutoff, row["threshold"])
            close(min(omissions.values()), row["loo_min"])
            close(max(omissions.values()), row["loo_max"])
            assert expected["stable_above"] + expected["sensitive"] + expected["stable_below"] == 199
            for omission in marker["leave_one_out"]:
                close(omission["mean"] + k * omission["sd"], omissions[omission["omitted_id"]])
            verified_scenarios += 1
    print(json.dumps({"status": "PASS", "independent_implementation": "Python csv + statistics, no pandas or numpy",
        "exact_source_values_checked": 976, "gene_threshold_scenarios": verified_scenarios,
        "leave_one_control_out_thresholds_checked": verified_scenarios * 45,
        "groups_disjoint": all(not set(cases[g]) & set(controls[g]) for g in cases)}, indent=2))

if __name__ == "__main__":
    main()
