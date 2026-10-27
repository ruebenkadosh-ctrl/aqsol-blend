from pathlib import Path
import sys
import numpy as np,joblib
from rdkit import Chem
from rdkit.Chem import rdReducedGraphs
from lightgbm import Booster
from catboost import CatBoostRegressor
from features import featurize

def predict(smiles,model_dir=None,cache_dir=None):
 root=Path(model_dir) if model_dir else Path(__file__).resolve().parent
 sys.path.insert(0,str(root/'mordred_dependencies'))
 from mordred import Calculator,descriptors
 values=list(smiles);x,meta=featurize(values,Path(cache_dir) if cache_dir else root/'feature_cache')
 if not values:return {'per_seed':np.empty((5,0)),'mean':np.empty(0)}
 mols=[Chem.MolFromSmiles(s) for s in values];erg=np.asarray([rdReducedGraphs.GetErGFingerprint(m,fuzzIncrement=.3,minPath=1,maxPath=15) for m in mols],dtype=np.float32);aug=np.concatenate([x,erg],axis=1)
 calc=Calculator([d for d in Calculator(descriptors,ignore_3D=True).descriptors if str(d)!='FilterItLogS'],ignore_3D=True)
 rows=[]
 for mol in mols:
  row=[]
  for v in calc(mol):
   try:v=float(v)
   except (TypeError,ValueError):v=float('nan')
   row.append(v if np.isfinite(v) and abs(v)<=np.finfo(np.float32).max else 0.)
  rows.append(row)
 mx=np.asarray(rows,dtype=np.float32);predictions=[]
 for seed in range(1,6):
  sv=joblib.load(root/f'seed{seed}_svr.joblib');sp=sv['model'].predict(np.clip(sv['scaler'].transform(x[:,sv['descriptor_start']:].astype(np.float64)),-10,10))
  lgb=Booster(model_file=str(root/f'seed{seed}_lgb.model'));cat=CatBoostRegressor();cat.load_model(str(root/f'seed{seed}_cat.model'))
  predictions.append(.675*lgb.predict(aug)+.225*cat.predict(mx)+.1*sp)
 result=np.asarray(predictions);assert np.isfinite(result).all();return {'per_seed':result,'mean':result.mean(axis=0)}
