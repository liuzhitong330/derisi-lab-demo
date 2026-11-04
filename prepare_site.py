#!/usr/bin/env python3
"""Create minimal browser assets and two static figures from verified outputs."""
import argparse
import html
import json
import shutil
from pathlib import Path

def svg_begin(title, description, width=680, height=320):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        f'<title>{html.escape(title)}</title><desc>{html.escape(description)}</desc>',
        '<style>text{font-family:Arial,sans-serif;font-size:12px;fill:#222}.heading{font-family:Georgia,serif;font-size:19px}.note{font-size:11px;fill:#555}</style>',
        '<rect width="100%" height="100%" fill="#fff"/>',
        f'<text class="heading" x="24" y="30">{html.escape(title)}</text>']

def write_figures(data, out):
    markers = data["genes"]
    lines = svg_begin("Which calls depend on one control?", "Four antigens, 199 MIS-C specimens each; mean plus three population control standard deviations; leave-one-control-out categories.")
    colors = ["#1f7a8c", "#fff", "#bbb"]
    for i, marker in enumerate(markers):
        row = marker["settings"]["3"]
        counts = [row["stable_above"], row["sensitive"], row["stable_below"]]
        y = 80 + i * 45
        lines.append(f'<text x="24" y="{y+15}">{marker["gene"]}</text>')
        x = 100.0
        for count, color in zip(counts, colors):
            w = count / 199 * 400
            lines.append(f'<rect x="{x}" y="{y}" width="{w}" height="22" fill="{color}" stroke="#555" stroke-width="0.5"/>')
            x += w
        lines.append(f'<text x="518" y="{y+15}">{counts[0]} / {counts[1]} / {counts[2]}</text>')
    lines.extend(['<text class="note" x="24" y="275">Teal: above all 45 references · white: changes · gray: above none</text>',
        '<text class="note" x="24" y="295">Descriptive control-reference categories, not confirmed binding or clinical labels.</text>', '</svg>'])
    (out / "figures" / "control-reference-stability.svg").write_text("\n".join(lines) + "\n")
    lines = svg_begin("How far does the reference move?", "All-control cutoff and leave-one-control-out range for four antigens at k equals three; processed gene-level fold change over mock-IP, linear axis.")
    left, right, xmax = 110, 620, 4.0
    position = lambda value: left + value / xmax * (right-left)
    for tick in range(5):
        x = position(tick)
        lines.append(f'<line x1="{x}" y1="60" x2="{x}" y2="245" stroke="#e5e5e5"/>')
        lines.append(f'<text x="{x}" y="263" text-anchor="middle">{tick}</text>')
    for i, marker in enumerate(markers):
        row, y = marker["settings"]["3"], 82 + i * 45
        lines.append(f'<text x="24" y="{y+4}">{marker["gene"]}</text>')
        lines.append(f'<line x1="{position(row["loo_min"])}" y1="{y}" x2="{position(row["loo_max"])}" y2="{y}" stroke="#aaa" stroke-width="8"/>')
        lines.append(f'<circle cx="{position(row["threshold"])}" cy="{y}" r="4" fill="#1f7a8c"/>')
        lines.append(f'<text class="note" x="{position(row["loo_max"])+10}" y="{y+4}">{row["loo_min"]:.3f}–{row["loo_max"]:.3f}</text>')
    lines.extend(['<text class="note" x="24" y="285">Point: all-control cutoff · line: omitted-control range, NOT a confidence interval.</text>',
        '<text class="note" x="24" y="305">Gene-level FC/mock-IP. A cutoff below 1 does not establish binding or enrichment.</text>', '</svg>'])
    (out / "figures" / "cutoff-perturbation.svg").write_text("\n".join(lines) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--derived", type=Path, default=Path("derived"))
    parser.add_argument("--out", type=Path, default=Path("."))
    args = parser.parse_args()
    data = json.loads((args.derived / "demo-data.json").read_text())
    assert data["verification"]["full_matrix_z_parity"]
    assert len(data["genes"]) == 4
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(exist_ok=True)
    (args.out / "results.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    (args.out / "data.js").write_text("window.UNIVERSITY_LAB_DEMO_DATA = " + json.dumps(data, separators=(",", ":"), allow_nan=False) + ";\n")
    shutil.copyfile(args.derived / "sample_reference_audit.csv", args.out / "specimen_calls.csv")
    shutil.copyfile(args.derived / "threshold_summary.csv", args.out / "threshold_summary.csv")
    shutil.copyfile(args.derived / "verification.json", args.out / "verification.json")
    write_figures(data, args.out)
    print("Prepared results.json, data.js, CSVs, verification and two SVG figures.")

if __name__ == "__main__":
    main()
