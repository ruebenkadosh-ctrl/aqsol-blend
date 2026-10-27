"""Stage an isolated workspace for the archived external workflow; never score."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--train-val', type=Path, required=True)
    parser.add_argument('--reference-test', type=Path, required=True,
                        help='AqSolDB reference CSV: only Drug is copied, never test labels')
    parser.add_argument('--bigsoldb', type=Path, required=True)
    parser.add_argument('--current-package', type=Path, required=True)
    parser.add_argument('--previous-package', type=Path, required=True)
    args = parser.parse_args()
    import pandas as pd
    archive = Path(__file__).resolve().parents[1] / 'analysis/external'
    expected = json.loads((archive / 'feasibility.json').read_text())['dataset_sha256']
    if hashlib.sha256(args.bigsoldb.read_bytes()).hexdigest() != expected:
        raise ValueError('BigSolDB input differs from the recorded source')
    # Validate both frozen packages before staging anything.
    for package in [args.current_package, args.previous_package]:
        for name, digest in json.loads((package / 'manifest.json').read_text()).items():
            if hashlib.sha256((package / name).read_bytes()).hexdigest() != digest:
                raise ValueError(f'Package integrity failure: {name}')
    args.output.mkdir(parents=True, exist_ok=False)
    data = args.output / 'work/tdc/data/admet_group/solubility_aqsoldb'
    data.mkdir(parents=True)
    shutil.copyfile(args.train_val, data / 'train_val.csv')
    pd.read_csv(args.reference_test, usecols=['Drug']).to_csv(data / 'test.csv', index=False)
    analysis = args.output / 'work/aqsol_external'
    analysis.mkdir(parents=True)
    for name in ['feasibility.py', 'refine_cohort.py', 'freeze_cohort.py',
                 'run_external.py', 'source_schema.json', 'external_protocol.json']:
        shutil.copyfile(archive / name, analysis / name)
    shutil.copyfile(args.bigsoldb, analysis / 'BigSolDBv2.2.csv')
    outputs = args.output / 'outputs'
    outputs.mkdir()
    for source, name in [(args.current_package, 'solubility_mordred_candidate'),
                         (args.previous_package, 'solubility_erg_svr_candidate')]:
        shutil.copytree(source, outputs / name, ignore=shutil.ignore_patterns('__pycache__', 'feature_cache'))
    (args.output / 'README.txt').write_text(
        'Independent reproduction workspace, not the original experiment.\n'
        'From this directory run, in order:\n'
        'python work/aqsol_external/feasibility.py\n'
        'python work/aqsol_external/refine_cohort.py\n'
        'python work/aqsol_external/freeze_cohort.py --eligible work/aqsol_external/eligible_metadata.csv '
        '--protocol work/aqsol_external/external_protocol.json '
        '--output work/aqsol_external/cohort_metadata_frozen.csv\n'
        'Only after metadata verification, the optional evaluation command is:\n'
        'python work/aqsol_external/run_external.py\n'
        'It trains the fixed baseline and scores the frozen models once. '
        'Results belong to this reproduction; they must not replace the original ledgers.\n')
    print('Workspace staged. No training, inference, or scoring performed.')


if __name__ == '__main__':
    main()
