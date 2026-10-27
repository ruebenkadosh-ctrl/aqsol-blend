"""Validate report consistency and Python syntax without training or rescoring."""
from pathlib import Path
import ast,hashlib,json
root=Path(__file__).resolve().parents[1]
for p in list((root/'src').glob('*.py'))+list((root/'scripts').glob('*.py'))+list((root/'analysis').rglob('*.py')):ast.parse(p.read_text(),filename=str(p))
result=root/'results/external_results.json';ledger=json.loads((root/'results/external_complete.json').read_text());assert hashlib.sha256(result.read_bytes()).hexdigest()==ledger['results_sha256']
a=json.loads(result.read_text());assert a['n_molecules']==360
b=json.loads((root/'results/official_results.json').read_text());assert b['solubility_aqsoldb']==[.736,.014]
c=json.loads((root/'config.json').read_text());assert abs(sum(c['weights'])-1)<1e-12
archive=root/'analysis/external'
for name,digest in json.loads((archive/'archive_manifest.json').read_text()).items():
 assert hashlib.sha256((archive/name).read_bytes()).hexdigest()==digest, name
print('PASS: source syntax, result hash, configuration, and archived source integrity. No model fit or score recomputation.')
