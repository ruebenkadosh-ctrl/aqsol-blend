"""Package trusted frozen models without datasets, validation labels, or caches."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    names = [f'seed{seed}_{branch}.{suffix}' for seed in range(1, 6)
             for branch, suffix in [('lgb', 'model'), ('cat', 'model'), ('svr', 'joblib')]]
    source_manifest = json.loads((args.source / 'manifest.json').read_text())
    manifest = {}
    for name in names:
        digest = sha256(args.source / name)
        if source_manifest.get(name) != digest:
            raise ValueError(f'Frozen artifact hash mismatch: {name}')
        manifest[name] = digest
    with zipfile.ZipFile(args.output, mode='x', compression=zipfile.ZIP_DEFLATED) as archive:
        for name in names:
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 26, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, (args.source / name).read_bytes())
        archive.writestr(zipfile.ZipInfo('manifest.json', (2026, 10, 26, 0, 0, 0)),
                         json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    with zipfile.ZipFile(args.output) as archive:
        for name, digest in manifest.items():
            if hashlib.sha256(archive.read(name)).hexdigest() != digest:
                raise ValueError(f'Archive verification failed: {name}')
    print(json.dumps({'archive': args.output.name, 'sha256': sha256(args.output),
                      'model_files': len(names), 'bytes': args.output.stat().st_size}, indent=2))


if __name__ == '__main__':
    main()
