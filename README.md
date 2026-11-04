# Which antibody-screen calls survive a different control reference?

[Open the interactive demo](https://liuzhitong330.github.io/derisi-lab-demo/). Prepared by Cathy Liu, November 2026, with AI assistance. Personal research exploration; no lab affiliation or endorsement.

## Result in brief

At a descriptive mean + 3 population-SD cutoff, SNX8 has 75 of 199 MIS-C specimens above the all-control boundary. Of all 199 specimens, 74 remain above every leave-one-control-out boundary, 12 switch category across those boundaries, and 113 remain below. The 12 are **not** 12 of the 75: sensitivity can move a call in either direction. Omitting the most influential control lowers the SNX8 threshold by 10.5%, but is not evidence that this control should be excluded.

| Antigen | Above all-control cutoff | Above every LOO cutoff | Reference-sensitive | Below every LOO cutoff |
|---|---:|---:|---:|---:|
| SNX8 | 75 | 74 | 12 | 113 |
| ERFL | 102 | 100 | 3 | 96 |
| KDELR1 | 70 | 70 | 4 | 125 |
| UBE3A | 42 | 42 | 10 | 147 |

These are **screening-reference categories, not confirmed binding**. A sensible follow-up is to test a balanced, blinded panel spanning stable-above, reference-sensitive and stable-below specimens on an orthogonal full-length protein assay, after checking treatment metadata and material availability. The original paper already includes orthogonal validation; this tool makes a narrower specimen-selection question inspectable.

Two static summaries: [call stability](figures/control-reference-stability.svg) and [cutoff movement](figures/cutoff-perturbation.svg).

## Reproduce

This is a focused reanalysis of processed measurements from Bodansky et al., *Nature* (2024), DOI [10.1038/s41586-024-07722-4](https://doi.org/10.1038/s41586-024-07722-4). It is not a reproduction of the paper's full feature-selection/classifier pipeline and not an independent validation cohort.

## Inputs

Download these two files from [Dryad DOI 10.7272/Q6SJ1HVH](https://datadryad.org/dataset/doi:10.7272/Q6SJ1HVH) into `inputs/` using the ordinary browser download links:

- `Human_Proteome_MIS-C_patients.csv`
- `Human_Proteome_at-risk_controls.csv`

The API metadata and SHA256 checksums are in `provenance.json`. API download requests currently require authentication and direct scripted file-stream requests returned an error; browser downloads worked. No account or private data is required for the ordinary public file links. The data are CC0. Do not substitute another version without reviewing its changes and updating the manifest deliberately.

Rows are human genes, columns are deidentified study sample IDs, and values are processed gene-level fold changes over mock-immunoprecipitation. Per the paper's Methods, peptide reads are summed by gene, a 0.5-read pseudocount is added, values are normalized to read percentages, then divided by the mean corresponding bead-only mock-IP read percentage. They are **not sums of peptide-level fold changes, raw counts, RPK values, antibody titres, or expression measurements**. The 199 MIS-C samples and 45 at-risk samples contain 19,997 matching gene rows. These children were previously infected; the controls are not pre-COVID healthy blood donors.

The focus is fixed in advance from the published work: SNX8, ERFL, KDELR1, and UBE3A as a control-reactivity comparator. The deposited row `LOC390937` is mapped explicitly to ERFL using [NCBI Gene 390937](https://www.ncbi.nlm.nih.gov/gene/390937); it is not ERF.

## Run

Use Python 3.11 or newer. The analysis was checked with the pinned NumPy and pandas versions below; only standard-library modules are used by the independent verifier.

```sh
python -m pip install -r requirements.txt
git clone https://github.com/h-s-miller/phagepy.git phagepy
git -C phagepy checkout d600a664c69541f8b04fe47fdce166fe52c16f6f
python analyze.py --inputs inputs --phagepy phagepy --out derived
python check_results.py --inputs inputs --derived derived
python prepare_site.py --derived derived --out .
```

The upstream checkout's `analysis.py` was read before execution. The wrapper verifies both the commit and that module's SHA256, then executes its unchanged `generate_SS_hits` function. It uses a deliberately small AnnData-like numeric adapter exposing only the objects used by that function. It does not invoke raw-count loading or renormalize the processed fold changes.

The phagepy repository is associated with the lab's published PhIP-seq methods, including [JCI Insight 174976](https://insight.jci.org/articles/view/174976). It is **not claimed to be the exact code used for the MIS-C data deposit**. No explicit software license was found in the pinned checkout, so its code is not redistributed here; it is fetched separately. This project's wrapper and verifier are new code.

## Decisions and outputs

The screening boundary is the at-risk control arithmetic mean plus `k` population standard deviations (`ddof=0`, matching the inspected function). A sample must be strictly greater than the boundary to be called “above reference.” Default `k=3`; the interface also compares 2, 4, and 5. These are descriptive sensitivity settings, not optimized or clinical thresholds.

For each antigen and setting, each of the 45 controls is omitted once and the boundary is recomputed. Case observations are above all 45 boundaries, above some, or above none. Each control is separately assessed against the other 44 controls. The leave-one-out band is a deterministic perturbation range—not a confidence interval or probability.

The paper already describes evaluating control z scores against the remaining controls. Our extension is the systematic control-omission sensitivity of **case** calls plus the interactive selection view; it is not the invention of held-out-control scoring.

Outputs:

- `derived/demo-data.json`: selected sample values, control-omission summaries, audit metadata and numerical verification.
- `derived/threshold_summary.csv`: 16 antigen/setting scenarios.
- `derived/sample_reference_audit.csv`: 199 case observations × 4 antigens × 4 settings, with reproducible reference categories.
- `derived/verification.json`: upstream-function parity checks and independent spot checks.
- `results.json`, `data.js`, and `specimen_calls.csv`: minimal browser/shareable assets generated from those outputs.
- `figures/control-reference-stability.svg` and `figures/cutoff-perturbation.svg`: two static summaries at the default setting.

`check_results.py` independently reads the original CSVs and uses Python's `statistics` module, not NumPy or pandas, to verify all 976 selected input values, the 16 scenarios and 720 leave-one-out boundaries.

## What the analysis can support

It can help propose a balanced **future assay challenge set**: robust-above-reference samples, boundary-sensitive samples, below-reference comparisons and high-control challenges. This is more informative for assay-method evaluation than only taking the largest case signals. Any shortlist requires lab confirmation of material availability, independent aliquots, treatment and collection-time metadata before experimentation. It is not a claim that these samples were never tested in the original study.

The paper already contains orthogonal validation of its main findings. A future assay comparison would need a predeclared panel, blinded sample order, blank and positive controls, an irrelevant-antigen control, dilution checks, and a genuinely independent batch or cohort. Discordance can reflect assay format or sample effects and does not itself refute a mechanism.

## Limits

- No individual IVIG status, collection timing, raw reads, technical replicates, aliquot inventory or linked orthogonal-assay results accompany these CSVs. Treatment-adjusted estimates are therefore not made.
- Mean/SD reference boundaries are sensitive to skewed control values. The audit displays that sensitivity; it does not establish a correct biological threshold.
- The KDELR1 default boundary is below 1 fold change over mock-IP. An above-at-risk-reference call must not be described as confirmed binding or enrichment relative to mock-IP. Do not introduce a universal FC>1 filter as a validated cutoff.
- Leaving out controls does not create a new cohort. The targets were chosen using the same published discovery cohort, so apparent separation is not out-of-sample diagnostic performance.
- High- or low-valued controls are not presumed bad samples. Omission is a diagnostic stress test, not an exclusion recommendation.
- The viral CSVs contain 181 cases and 45 controls, despite the README's generic 199-case description. They were checked for provenance but are not used here; the two libraries have different normalization denominators.

This artifact was developed with AI coding assistance. It demonstrates a transparent research workflow and experiment-planning approach, not independent completion of PhIP-seq, primary-cell culture, cloning, or library preparation.
