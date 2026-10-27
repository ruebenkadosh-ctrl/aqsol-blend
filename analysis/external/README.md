# External evaluation source archive

These are unmodified source scripts and provenance records from the completed evaluation. `archive_manifest.json` records their SHA-256 digests. They preserve the original computation rather than representing a portable, independently verified reproduction command.

## Original execution layout

The original project placed these scripts in `work/aqsol_external/` and ran cohort preparation from its root. Inputs were:

- `work/tdc/data/admet_group/solubility_aqsoldb/train_val.csv`;
- `work/tdc/data/admet_group/solubility_aqsoldb/test.csv` (only the `Drug` structure column for overlap removal);
- `outputs/solubility_mordred_candidate/`, the current frozen package;
- `outputs/solubility_erg_svr_candidate/`, the preceding frozen package;
- BigSolDB v2.2, obtained from the pinned source in `source_schema.json` and checked against `feasibility.json`.

Both model packages contain their original prediction code and file manifests. They are not yet distributed here. Raw datasets and molecule-level outputs are also not included in this archive.

## Sequence and reconstruction check

1. `feasibility.py` reads only molecular structures and metadata to construct the initial water/temperature/identity filter.
2. `refine_cohort.py` conservatively excludes structures matching the heavy-element composition of unresolved reference structures and computes nearest-training similarity.
3. `freeze_cohort.py` reconstructs the originally separate metadata-only freeze operation: exclude fingerprint Tanimoto 1.0, preserve row order and CSV columns, and verify the protocol digest before writing. It was checked against the original intermediate metadata and recreated the exact frozen digest (421 rows, 361 parents). Unlike the three archived scripts, this is a new reconstruction utility, not original execution evidence. Use `--eligible eligible_metadata.csv --protocol external_protocol.json --output cohort_metadata_frozen.csv` in a separate reproduction directory. Raw inputs and frozen model distribution remain pending.
4. `run_external.py` verifies the frozen cohort and model manifests, trains the fixed RF baseline, saves predictions, and only then loads external labels for scoring.

The original start and completion ledgers are preserved. The runner refuses to overwrite its start marker. Do not delete that marker to rescore this completed experiment. An independent reproduction must use a separate output location and distinguish itself from the original evaluation.

The protocol's status text records its state when frozen. The completed evaluation and source review document the subsequent execution; the protocol has not been rewritten after observing results.

## Portable staging

`scripts/prepare_external_workspace.py` accepts explicit input and model-package paths and creates the original directory layout in a new destination. It validates input/model hashes, copies only `Drug` from the AqSolDB test reference, and writes ordered reproduction commands to the destination README. It does not train or score. Run its `--help` for required arguments. Public input/model distribution remains separate.
