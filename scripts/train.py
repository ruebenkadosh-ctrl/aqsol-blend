"""Fixed development training recipe. Does not score official test labels."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='6'
import argparse,json,sys
from pathlib import Path
import numpy as np,joblib
from rdkit import Chem
from rdkit.Chem import rdReducedGraphs
from mordred import Calculator,descriptors
from lightgbm import LGBMRegressor,early_stopping
from catboost import CatBoostRegressor
from sklearn.preprocessing import RobustScaler
from sklearn.svm import SVR
from tdc.benchmark_group import admet_group
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
from features import featurize

def main():
 p=argparse.ArgumentParser();p.add_argument('--data-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,choices=range(1,6),required=True);a=p.parse_args();cfg=json.loads((root/'config.json').read_text());out=a.output/f'seed{a.seed}';out.mkdir(parents=True,exist_ok=False)
 tr,va=admet_group(path=str(a.data_root)).get_train_valid_split(seed=a.seed,benchmark='solubility_aqsoldb',split_type='default')
 smiles=tr.Drug.tolist()+va.Drug.tolist();n=len(tr);x,meta=featurize(smiles,out/'cache');assert x.shape[1]==2937,'Feature schema mismatch'
 mols=[Chem.MolFromSmiles(s) for s in smiles];erg=np.asarray([rdReducedGraphs.GetErGFingerprint(m,fuzzIncrement=.3,minPath=1,maxPath=15) for m in mols],dtype=np.float32);aug=np.concatenate([x,erg],axis=1)
 calc=Calculator([d for d in Calculator(descriptors,ignore_3D=True).descriptors if str(d)!='FilterItLogS'],ignore_3D=True);rows=[]
 for mol in mols:
  row=[]
  for value in calc(mol):
   try:value=float(value)
   except (TypeError,ValueError):value=float('nan')
   row.append(value if np.isfinite(value) and abs(value)<=np.finfo(np.float32).max else 0.)
  rows.append(row)
 mx=np.asarray(rows,dtype=np.float32);assert mx.shape[1]==1612
 l=LGBMRegressor(**cfg['lightgbm'],random_state=a.seed);l.fit(aug[:n],tr.Y,eval_set=[(aug[n:],va.Y)],callbacks=[early_stopping(150,verbose=False)]);l.booster_.save_model(str(out/f'seed{a.seed}_lgb.model'))
 c=CatBoostRegressor(**cfg['catboost'],grow_policy='SymmetricTree',random_seed=a.seed);c.fit(mx[:n],tr.Y,eval_set=(mx[n:],va.Y),early_stopping_rounds=150,verbose=False);c.save_model(str(out/f'seed{a.seed}_cat.model'))
 start=meta['groups']['descriptors'][0];xd=x[:,start:].astype(np.float64);sc=RobustScaler(quantile_range=(10,90)).fit(xd[:n]);s=SVR(**cfg['svr']);s.fit(np.clip(sc.transform(xd[:n]),-10,10),tr.Y);joblib.dump({'scaler':sc,'model':s,'descriptor_start':start,'clip':10},out/f'seed{a.seed}_svr.joblib')
 pred=.675*l.predict(aug[n:])+.225*c.predict(mx[n:])+.1*s.predict(np.clip(sc.transform(xd[n:]),-10,10));np.savez_compressed(out/'validation.npz',smiles=va.Drug.to_numpy(dtype=str),truth=va.Y.to_numpy(),prediction=pred)
 (out/'summary.json').write_text(json.dumps({'seed':a.seed,'validation_mae':float(np.mean(abs(va.Y.to_numpy()-pred))),'disclosure':'Outer-validation early stopping, as in original recipe. Not nested confirmation.'},indent=2))
if __name__=='__main__':main()
