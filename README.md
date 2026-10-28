# AqSol-Blend

**A CPU-Trained Descriptor Ensemble for Aqueous Solubility Prediction**  
Rueben Kadosh

[Paper](docs/manuscript.pdf) · [Trained models](https://github.com/ruebenkadosh-ctrl/aqsol-blend/releases/tag/v0.1.0) · [Reproducibility](REPRODUCIBILITY.md)

AqSol-Blend combines LightGBM, CatBoost and support-vector regression to predict aqueous solubility from molecular structure. This repository accompanies the manuscript and provides training and inference scripts, frozen model files, and evaluation results. The experiments examine how molecular representations affect an ensemble of conventional regressors, with evaluation on TDC AqSolDB and an overlap-screened BigSolDB cohort.

## Results

| Evaluation | Model | MAE | RMSE |
| --- | --- | ---: | ---: |
| TDC AqSolDB | AqSol-Blend | 0.736 ± 0.014 | — |
| BigSolDB external cohort (360 compounds) | AqSol-Blend | 1.008 | 1.318 |
| Same external cohort | Previous ensemble | 1.017 | 1.335 |
| Same external cohort | Fixed RDKit random forest | 1.065 | 1.415 |

TDC values are the mean and standard deviation of five per-seed MAEs. External predictions use the average of the five ensembles. On the external cohort, the paired MAE improvement over random forest is 0.058 (95% bootstrap interval 0.022–0.091); the improvement over the preceding ensemble is 0.009 (−0.001–0.019).

### Evaluation scope

Four sequential TDC test evaluations occurred during development, so the benchmark result reflects adaptive test exposure. It is a locally computed result, not an accepted leaderboard ranking. The external protocol and predictions were fixed before scoring. BigSolDB measurements differ in concentration convention and available measurement metadata; the manuscript describes these limitations and the overlap-removal procedure.

## Model

Each seed blends:

- 67.5% LightGBM using fingerprints, RDKit descriptors and ErG features;
- 22.5% CatBoost using 1,612 Mordred descriptors;
- 10% RBF SVR using robust-scaled RDKit descriptors.

The mixture weights were fixed before external evaluation. The manuscript reports development comparisons, branch ablations and weight sensitivity.

## Environment

Use Python 3.11 and install `pip install -r requirements.txt` in a virtual environment. Native LightGBM may require OpenMP on macOS. The pinned requirements installed successfully in a fresh Python 3.11 macOS ARM64 environment, and `pip check` reported no conflicts. The resolved dependency list is in `requirements-resolved-macos-arm64.txt`. All five seeds were retrained in that environment; their saved validation predictions matched the originals exactly. RDKit version changes can alter descriptor values and ordering.

## Inference

Download the fifteen frozen model files from [release v0.1.0](https://github.com/ruebenkadosh-ctrl/aqsol-blend/releases/tag/v0.1.0). Extract the archive into one directory and use that directory as `--models`:

```sh
python scripts/predict_smiles.py --models /path/to/frozen/package --smiles CCO --output predictions.json
```

The output contains five per-seed predictions and their mean. Only load trusted joblib artifacts.

## Training

```sh
python scripts/train.py --data-root data --output training_runs --seed 1
```

Repeat for seeds 1–5. Each run writes a separate `training_runs/seedN` directory and refuses to overwrite it. The inference loader expects the fifteen model files together in one model directory; assemble them explicitly after training. The training script does not score official test labels. Fresh-environment retraining on the same host reproduced all five sets of validation predictions exactly; see `results/training_reproduction.json`. Descriptor timeouts may depend on hardware.

## Evaluation artifacts

- `analysis/external/`: original evaluation scripts, provenance, and documented portability gaps.
- `results/official_results.json`: saved TDC evaluator output.
- `results/external_protocol.json`: fixed external evaluation protocol.
- `results/external_results.json`: external metrics and paired comparisons.
- `results/external_complete.json`: external result digest.
- `results/branch_analysis.json`: development branch ablations.
- `results/weight_sensitivity.json`: descriptive sensitivity on reused development splits; weights were not changed.
- `results/development_checks.json`: development sanity checks.
- `results/model_archive.json`: model archive checksum and size.
- `results/feature_schema.json`: ordered descriptors and feature settings.
- `results/MODEL_INVENTORY.json`: trees, leaves, support vectors and artifact sizes.

Run `python scripts/check_release.py` for report consistency and syntax checks. This does not retrain models or recompute scores.

## Data availability

Obtain AqSolDB through TDC and BigSolDB v2.2 from [Zenodo](https://zenodo.org/records/22648301). Raw datasets are not distributed in this repository. BigSolDB is distributed under CC BY 4.0; other datasets and dependencies retain their respective licenses.

## License and manuscript

Original project code is released under the [MIT License](LICENSE). The [manuscript](docs/manuscript.pdf) is available under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). It has not yet been announced on arXiv. For author and citation metadata, see [CITATION.cff](CITATION.cff).
