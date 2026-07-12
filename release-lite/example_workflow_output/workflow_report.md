# Tumor-vs-normal workflow report

Input: `example_input.csv`

## QC

- Status: **PASS**
- Model genes matched: 2000/2000 (100.0%)
- Expression median / p99 / max: 3.6067 / 9.0798 / 13.8937
- |z| > 6 fraction: 0.001100
- Cohort gene-mean |z| p99: 1.9613

No QC warnings or errors.

## Scores

- Samples: 5
- Tumor calls: 3
- Normal calls: 2
- Tumor probability median / p90 / max: 0.9900 / 0.9999 / 0.9999

## Calibration

- Labeled samples: 5 (3 tumor / 2 normal)
- AUC: 1.0000
- Evaluation: **apparent/resubstitution** (threshold and metrics were estimated on these same labeled samples; this is not independent validation)
- Recommended threshold: 0.989957 (youden_j)
- Recommended accuracy / recall / specificity: 1.0000 / 1.0000 / 1.0000

- WARNING: Calibration metrics are unstable because at least one class has fewer than 10 labeled samples (3 tumor, 2 normal).

## Explanations

- Explanation rows: 30
- Positive rows push the LR logit toward tumor; negative rows push it toward normal.

## Output files

| file | path |
| --- | --- |
| qc_json | qc.json |
| scores_csv | scores.csv |
| report_md | workflow_report.md |
| manifest_json | manifest.json |
| thresholds_csv | thresholds.csv |
| calibration_json | calibration.json |
| explanations_csv | explanations.csv |
