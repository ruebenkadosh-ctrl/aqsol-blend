import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import json,collections
import pandas as pd,numpy as np
from rdkit import Chem
from mordred import Calculator,descriptors
from mordred.error import MissingValueBase,Timeout
r=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
calc=Calculator([d for d in Calculator(descriptors,ignore_3D=True).descriptors if str(d)!='FilterItLogS'],ignore_3D=True)
sets={'test':pd.read_csv(r/'work/tdc/data/admet_group/solubility_aqsoldb/test.csv',usecols=['Drug']).Drug.tolist(),'external':sorted(pd.read_csv(r/'work/aqsol_external/cohort_metadata_frozen.csv',usecols=['SMILES_Solute']).SMILES_Solute.unique())}
for name,smiles in sets.items():
 path=out/f'{name}_descriptor_audit.jsonl'
 with path.open('x') as f:
  for i,s in enumerate(smiles):
   failures=collections.Counter();bad=[]
   for d,v in zip(calc.descriptors,calc(Chem.MolFromSmiles(s))):
    reason=None
    if isinstance(v,MissingValueBase):reason='timeout' if isinstance(v.error,Timeout) else type(v.error).__name__
    else:
     try:
      z=float(v)
      if not np.isfinite(z):reason='nonfinite'
      elif abs(z)>np.finfo(np.float32).max:reason='float32_overflow'
     except (TypeError,ValueError):reason='conversion_error'
    if reason:failures[reason]+=1;bad.append(str(d))
   f.write(json.dumps({'row':i,'failures':dict(failures),'affected_descriptors':bad})+'\n');f.flush()
   if (i+1)%100==0:print(name,i+1,'/',len(smiles),flush=True)
 rows=[json.loads(line) for line in path.read_text().splitlines()];counts=collections.Counter()
 for row in rows:counts.update(row['failures'])
 result={'scope':'Fresh label-free extraction audit; not historical timeout incidence','n_input_rows':len(rows),'descriptors_per_row':len(calc.descriptors),'molecules_with_any_failure':sum(bool(x['failures']) for x in rows),'molecules_with_timeout':sum(x['failures'].get('timeout',0)>0 for x in rows),'descriptor_failure_counts':dict(counts)}
 (out/f'{name}_descriptor_summary.json').write_text(json.dumps(result,indent=2));print(result,flush=True)
