"""Compare aligned development predictions; no official test or external labels."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', type=Path, required=True)
    p.add_argument('--retrained', type=Path, required=True)
    p.add_argument('--train-val', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    frame = pd.read_csv(a.train_val)
    results = []
    for seed in range(1, 6):
        with np.load(a.original / f'seed{seed}_validation.npz', allow_pickle=False) as original, \
             np.load(a.retrained / f'seed{seed}/validation.npz', allow_pickle=False) as new:
            expected = frame.iloc[original['feature_row']].Drug.to_numpy(dtype=str)
            if not np.array_equal(expected, new['smiles']):
                raise ValueError(f'Seed {seed}: molecule order mismatch')
            if not np.array_equal(original['truth'], new['truth']):
                raise ValueError(f'Seed {seed}: validation truth mismatch')
            delta = new['prediction'] - original['prediction']
            if not np.isfinite(delta).all():
                raise ValueError(f'Seed {seed}: nonfinite difference')
            results.append({'seed': seed, 'rows': len(delta),
                            'max_absolute_prediction_difference': float(np.max(np.abs(delta))),
                            'mean_absolute_prediction_difference': float(np.mean(np.abs(delta))),
                            'matches_at_1e_8': bool(np.allclose(new['prediction'], original['prediction'], atol=1e-8, rtol=0))})
    result = {'scope': 'Aligned saved development prediction comparison, not test scoring',
              'seeds': results, 'all_match_at_1e_8': all(x['matches_at_1e_8'] for x in results)}
    with a.output.open('x') as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
