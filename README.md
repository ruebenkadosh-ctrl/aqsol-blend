# AqSol-Blend

A CPU-trained ensemble for aqueous solubility prediction using LightGBM, CatBoost and support-vector regression. Research code release by **Rueben Kadosh**. This is not an accepted TDC submission.

## Results

| Evaluation | Model | MAE | RMSE |
| --- | --- | ---: | ---: |
| TDC AqSolDB | AqSol-Blend | 0.736 ± 0.014 | — |
| BigSolDB external cohort (360 compounds) | AqSol-Blend | 1.008 | 1.318 |
| Same external cohort | Previous ensemble | 1.017 | 1.335 |
| Same external cohort | Fixed RDKit random forest | 1.065 | 1.415 |

The TDC result is the mean and standard deviation of five per-seed MAEs. Four official evaluations occurred during adaptive development. It is not an independently confirmed first-place ranking.

External predictions average the five per-seed ensembles. The protocol was fixed and predictions committed before external labels were loaded for scoring. The paired bootstrap MAE improvement over random forest was 0.058 (95% interval 0.022–0.091); improvement over the previous ensemble was 0.009 (−0.001–0.019). The latter interval includes zero.

The external cohort uses BigSolDB v2.2 water measurements near 25 °C, with standardized molecular overlaps with AqSolDB removed. Its concentration convention, missing pH and solid-form information, and possible shared literature sources limit comparisons across datasets. It tests transfer to this cohort, not universal solubility generalization.

## Model

Each seed blends:

- 67.5% LightGBM using fingerprints, RDKit descriptors and ErG features;
- 22.5% CatBoost using 1,612 Mordred descriptors;
- 10% RBF SVR using robust-scaled RDKit descriptors.

The weights were inherited and fixed before external evaluation. These are established algorithms. The study documents their development, component contributions and external performance.

Read the [manuscript (PDF)](docs/manuscript.pdf). This manuscript has not yet been accepted or announced by arXiv. The manuscript is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the project code remains MIT-licensed. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) separates completed verification from work still pending.

## Environment

Use Python 3.11 and install `pip install -r requirements.txt` in a virtual environment. Native LightGBM may require OpenMP on macOS. The pinned requirements installed successfully in a fresh Python 3.11 macOS ARM64 environment, and `pip check` reported no conflicts. The resolved dependency list is in `requirements-resolved-macos-arm64.txt`. All five seeds were retrained in that environment; their saved validation predictions matched the originals exactly. RDKit version changes can alter descriptor values and ordering.

## Inference

A verified archive of the fifteen frozen model files has been prepared locally; the archive is available in [GitHub Releases](https://github.com/ruebenkadosh-ctrl/aqsol-blend/releases/tag/v0.1.0). Extract the archive into one directory and use that directory as `--models`:

```sh
python scripts/predict_smiles.py --models /path/to/frozen/package --smiles CCO --output predictions.json
```

The output contains five per-seed predictions and their mean. Only load trusted joblib artifacts.

## Training

```sh
python scripts/train.py --data-root data --output training_runs --seed 1
```

Repeat for seeds 1–5. Each run writes a separate `training_runs/seedN` directory and refuses to overwrite it. The inference loader expects the fifteen model files together in one model directory; assemble them explicitly after training. The training script does not score official test labels. Fresh-environment retraining on the same host reproduced all five sets of validation predictions exactly; see `results/training_reproduction.json`. Descriptor timeouts may depend on hardware.

## Recorded evidence

- `analysis/external/`: original evaluation scripts, provenance, and documented portability gaps.
- `results/official_results.json`: saved TDC evaluator output.
- `results/external_protocol.json`: fixed external evaluation protocol.
- `results/external_results.json`: external metrics and paired comparisons.
- `results/external_complete.json`: external result digest.
- `results/branch_analysis.json`: development branch ablations.
- `results/weight_sensitivity.json`: descriptive sensitivity on reused development splits; weights were not changed.
- `results/development_checks.json`: development sanity checks.
- `results/model_archive.json`: verified local model archive digest and size.
- `results/feature_schema.json`: ordered descriptors and feature settings.
- `results/MODEL_INVENTORY.json`: trees, leaves, support vectors and artifact sizes.

Run `python scripts/check_release.py` for report consistency and syntax checks. This does not retrain models or recompute scores.

## Release status and data availability

Author: Rueben Kadosh. Final manuscript review remains pending. Original project code is licensed under the MIT License; see LICENSE. This does not relicense third-party datasets, dependencies, or the manuscript. No credentials or raw benchmark data are included. Obtain AqSolDB through TDC and BigSolDB v2.2 from [its archived release](https://zenodo.org/records/22648301). BigSolDB is licensed CC BY 4.0; third-party software and other datasets retain their respective licenses.
