"""Reconstruct the metadata-only freeze step; verify before writing a new file."""
import argparse
import hashlib
import io
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eligible', required=True, type=Path)
    parser.add_argument('--protocol', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    frame = pd.read_csv(args.eligible)
    cohort = frame.loc[frame.max_train_tanimoto < 1.0].copy()
    protocol = json.loads(args.protocol.read_text())
    buffer = io.StringIO()
    cohort.to_csv(buffer, index=False)
    content = buffer.getvalue().encode('utf-8')
    digest = hashlib.sha256(content).hexdigest()
    if digest != protocol['cohort_sha256']:
        raise ValueError('Reconstructed metadata digest differs from the frozen protocol')
    if cohort.parent_identity.nunique() != protocol['n_parent']:
        raise ValueError('Parent count differs from the frozen protocol')
    with args.output.open('xb') as handle:
        handle.write(content)
    print(f'Verified {len(cohort)} rows; {protocol["n_parent"]} parents; SHA256 {digest}')


if __name__ == '__main__':
    main()
