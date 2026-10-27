"""Label-independent molecular features with row-order-preserving disk caching.

featurize(smiles, cache_dir) -> (X, metadata), where X is float32 and
metadata['groups'][name] is a JSON-compatible [start, stop] column slice.
Descriptors that fail or yield NaN/inf/out-of-float32-range values become 0;
metadata records the count. ECFP counts remain numeric counts (never uint8).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem, Descriptors, MACCSkeys

try:
    from rdkit.Avalon import pyAvalonTools
except ImportError:
    pyAvalonTools = None

FEATURE_VERSION = 'molecular-v1'
_DESCRIPTOR_LIST = tuple(Descriptors._descList)


def _spec():
    return {
        'version': FEATURE_VERSION,
        'rdkit_version': rdBase.rdkitVersion,
        'ecfp_radius': 2,
        'ecfp_columns': 2048,
        'avalon_columns': 512 if pyAvalonTools is not None else 0,
        'descriptors': [name for name, _ in _DESCRIPTOR_LIST],
        'nonfinite_descriptor_policy': 'replace failed, nonfinite, or float32-overflow values with zero',
    }


def _compute(smiles, spec):
    names = [f'ecfp_count_r2_{i}' for i in range(2048)]
    groups = {'ecfp': [0, len(names)]}
    start = len(names)
    names.extend(f'maccs_{i}' for i in range(167))
    groups['maccs'] = [start, len(names)]
    if pyAvalonTools is not None:
        start = len(names)
        names.extend(f'avalon_{i}' for i in range(512))
        groups['avalon'] = [start, len(names)]
    start = len(names)
    names.extend(f'rdkit_{name}' for name, _ in _DESCRIPTOR_LIST)
    groups['descriptors'] = [start, len(names)]
    groups['all'] = [0, len(names)]
    matrix = np.zeros((len(smiles), len(names)), dtype=np.float32)
    sanitized = 0
    limit = np.finfo(np.float32).max
    for row, value in enumerate(smiles):
        mol = Chem.MolFromSmiles(value)
        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError(f'Invalid or empty SMILES at index {row}: {value!r}')
        counts = AllChem.GetHashedMorganFingerprint(mol, radius=2, nBits=2048)
        for column, count in counts.GetNonzeroElements().items():
            matrix[row, column] = count
        a, b = groups['maccs']
        DataStructs.ConvertToNumpyArray(MACCSkeys.GenMACCSKeys(mol), matrix[row, a:b])
        if pyAvalonTools is not None:
            a, b = groups['avalon']
            DataStructs.ConvertToNumpyArray(pyAvalonTools.GetAvalonFP(mol, nBits=512), matrix[row, a:b])
        a, _ = groups['descriptors']
        for offset, (_, descriptor) in enumerate(_DESCRIPTOR_LIST):
            try:
                result = float(descriptor(mol))
            except Exception:
                result = float('nan')
            if not np.isfinite(result) or abs(result) > limit:
                sanitized += 1
                result = 0.0
            matrix[row, a + offset] = result
    return matrix, {
        **spec, 'feature_names': names, 'groups': groups,
        'n_rows': len(smiles), 'n_features': len(names),
        'sanitized_descriptor_values': sanitized,
    }


def featurize(smiles, cache_dir):
    """Compute/load features for an ordered iterable of SMILES strings.

    Keys hash exact input strings plus the full feature specification. Equivalent
    alternative SMILES therefore have separate cache entries, avoiding accidental
    row reordering or deduplication. Writes use atomic replacement; no labels or
    fitted statistics enter the features. Empty input yields a (0, n_features) array.
    """
    if isinstance(smiles, (str, bytes)):
        raise TypeError('smiles must be an iterable of strings, not one string')
    values = list(smiles)
    for i, value in enumerate(values):
        if not isinstance(value, str):
            raise ValueError(f'SMILES at index {i} must be a string')
    spec = _spec()
    encoded = json.dumps({'spec': spec, 'smiles': values}, ensure_ascii=True,
                         sort_keys=True, separators=(',', ':')).encode('utf-8')
    key = hashlib.sha256(encoded).hexdigest()
    directory = Path(cache_dir)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f'{key}.npz'
    if target.exists():
        with np.load(target, allow_pickle=False) as cached:
            matrix = cached['X']
            metadata = json.loads(str(cached['metadata'].item()))
        if (metadata.get('cache_key') != key or matrix.dtype != np.float32
                or matrix.shape != (len(values), metadata['n_features'])
                or not np.isfinite(matrix).all()):
            raise ValueError(f'Invalid feature cache: {target}')
        return matrix, metadata
    matrix, metadata = _compute(values, spec)
    metadata['cache_key'] = key
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(dir=directory, suffix='.npz', delete=False) as stream:
            temp_path = Path(stream.name)
            np.savez_compressed(stream, X=matrix, metadata=json.dumps(metadata, sort_keys=True))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_path, target)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
    return matrix, metadata
