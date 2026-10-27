# Reproducibility status

## Completed

- The frozen model package produced finite predictions in an inference smoke check.
- Fresh-environment ethanol inference reproduced all five original predictions within 1e-10 after the documented OpenMP configuration; see `results/clean_environment_smoke.json`. This single-molecule check does not verify full feature stability or retraining.
- The external evaluation completed under a fixed protocol, with predictions committed before label-based scoring.
- The saved external result digest matches its completion ledger.
- The manuscript and repository report the same rounded external MAE and RMSE values.
- Release Python files pass syntax parsing.
- The full metadata preparation workflow in a separate staged workspace reproduced the exact original cohort SHA-256 digest from source inputs (421 rows, 361 parents), without scoring models. See `results/cohort_reproduction.json`.
- Development ablations and weight sensitivity use saved, aligned predictions without new fitting.

All five seeds were retrained in the fresh environment on the same host. All 6,936 validation appearances matched the originals exactly after molecule/truth alignment. See `results/training_reproduction.json`. This is computational reproduction, not a new generalization evaluation or independent cross-platform replication.

## Still required for a complete public release

1. Complete runtime checks across the intended platforms. A fresh Python 3.11 macOS ARM64 install succeeded and `pip check` passed; the resolved dependency list is recorded. LightGBM initially failed to locate OpenMP on this host.
2. Independent replication on another host remains optional further evidence; the same-host fresh-environment five-seed retraining and aligned prediction comparison have passed.
3. Publish the verified model archive at a stable download location. Local packaging is complete: fifteen model files, each checked against its frozen source digest; archive digest is recorded in `results/model_archive.json`. Ordered descriptors are recorded in `results/feature_schema.json`.
4. Finish portability of the archived external workflow in `analysis/external/`, including input artifact distribution. The original scripts and provenance are now archived, but the complete end-to-end workflow is not yet portable.
5. Complete a final methodological and licensing review, author details, public links and the pending contribution/disclosure wording.

## Evaluation conventions

TDC results average five per-seed MAEs. External results score the mean prediction across five seeds. These quantities are different and must not be interchanged.

Development splits were reused and official test results were observed during development. The external comparison was frozen before its scoring, but endpoint and source differences remain. Further tuning on this external cohort would make it development data rather than a fresh evaluation.

## Artifact handling

Keep original models, prediction files and scoring ledgers immutable. Training creates separate per-seed directories. Never load untrusted joblib files. Do not commit credentials, local absolute paths, or raw data without checking redistribution rights.

## macOS OpenMP setup

The clean environment installed successfully but LightGBM could not initially find `libomp.dylib`. This host uses the OpenMP runtime bundled with the pinned scikit-learn wheel via `DYLD_LIBRARY_PATH`. For an active Python 3.11 virtual environment with that wheel layout:

```sh
export DYLD_LIBRARY_PATH="$VIRTUAL_ENV/lib/python3.11/site-packages/sklearn/.dylibs${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
```

This is a platform-specific runtime workaround, not an additional model dependency or a claim of Linux/Windows verification.

## Descriptor audit

A fresh label-free extraction on 1,997 test rows and 364 unique raw external SMILES recorded zero timeout-affected molecules. Other undefined descriptor cells numbered 383,761 and 78,735 respectively; all inputs had at least one such value. These use the frozen zero-replacement policy. Counts are preserved in `results/test_descriptor_summary.json` and `results/external_descriptor_summary.json`. They do not recover historical timeout incidence.
