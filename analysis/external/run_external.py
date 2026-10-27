import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='6'
from pathlib import Path
import json,hashlib,datetime,importlib.util,sys
import numpy as np,pandas as pd,joblib
from sklearn.ensemble import RandomForestRegressor
r=Path(__file__).resolve().parent;root=r.parents[1];h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,o):
 with p.open('x') as f:json.dump(o,f,indent=2,allow_nan=False)
def main():
 protocol=json.loads((r/'external_protocol.json').read_text());assert h(r/'cohort_metadata_frozen.csv')==protocol['cohort_sha256']
 write(r/'external_started.json',{'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':h(r/'external_protocol.json'),'runner_sha256':h(Path(__file__))})
 packages=[root/'outputs/solubility_mordred_candidate',root/'outputs/solubility_erg_svr_candidate']
 for p in packages:
  for n,s in json.loads((p/'manifest.json').read_text()).items():assert h(p/n)==s,(p,n)
 e=pd.read_csv(r/'cohort_metadata_frozen.csv');smiles=sorted(e.SMILES_Solute.unique());sys.path.insert(0,str(packages[0]));from features import featurize
 x,meta=featurize(smiles,r/'cache');dev=pd.read_csv(root/'work/tdc/data/admet_group/solubility_aqsoldb/train_val.csv');xt,mt=featurize(dev.Drug.tolist(),r/'cache');assert meta['feature_names']==mt['feature_names'];start=meta['groups']['descriptors'][0]
 cfg=protocol['baseline'];rf=RandomForestRegressor(**{k:cfg[k] for k in ['n_estimators','min_samples_leaf','max_features','random_state','n_jobs']});rf.fit(xt[:,start:],dev.Y);joblib.dump(rf,r/'rf.joblib');pred=pd.DataFrame({'smiles':smiles,'rf':rf.predict(x[:,start:])});print('RF baseline fitted',flush=True)
 for p,name in zip(packages,['current','previous']):
  spec=importlib.util.spec_from_file_location('inference_'+name,p/'predict.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);result=m.predict(smiles,model_dir=p,cache_dir=r/'cache');pred[name]=result['mean'];np.savez_compressed(r/(name+'_per_seed.npz'),smiles=np.array(smiles),predictions=result['per_seed']);print(name+' predictions saved',flush=True)
 assert np.isfinite(pred[['rf','current','previous']]).all().all();pred.to_csv(r/'external_predictions.csv',index=False);write(r/'predictions_committed.json',{'prediction_sha256':h(r/'external_predictions.csv')})
 # First target access, after all predictions have been saved.
 expected=json.loads((r/'feasibility.json').read_text())['dataset_sha256'];assert h(r/'BigSolDBv2.2.csv')==expected
 labels=pd.read_csv(r/'BigSolDBv2.2.csv',usecols=['SMILES_Solute','LogS(mol/L)']);assert (labels.iloc[e.source_row].SMILES_Solute.to_numpy()==e.SMILES_Solute.to_numpy()).all();e['truth']=labels.iloc[e.source_row]['LogS(mol/L)'].to_numpy();valid=np.isfinite(e.truth);missing=int((~valid).sum());e=e[valid]
 truth=e.groupby(['parent_identity','Source']).truth.median().groupby('parent_identity').median();ep=e[['parent_identity','SMILES_Solute']].drop_duplicates().merge(pred,left_on='SMILES_Solute',right_on='smiles');out=ep.groupby('parent_identity')[['current','previous','rf']].mean().join(truth).join(e.groupby('parent_identity').max_train_tanimoto.mean());out.to_csv(r/'external_molecule_results.csv')
 errors={k:np.abs(out[k].to_numpy()-out.truth.to_numpy()) for k in ['current','previous','rf']};metrics={k:{'mae':float(v.mean()),'rmse':float(np.sqrt(np.mean((out[k]-out.truth)**2)))} for k,v in errors.items()};rng=np.random.default_rng(20261026);idx=rng.integers(0,len(out),size=(2000,len(out)));paired={k:{'gain_current':float((errors[k]-errors['current']).mean()),'bootstrap95':np.quantile((errors[k]-errors['current'])[idx].mean(axis=1),[.025,.975]).tolist()} for k in ['previous','rf']};subsets={}
 for name,mask in [('below_0.3',out.max_train_tanimoto<.3),('0.3_to_0.6',(out.max_train_tanimoto>=.3)&(out.max_train_tanimoto<.6)),('at_least_0.6',out.max_train_tanimoto>=.6)]:
  subsets[name]={'n':int(mask.sum()),'mae':{k:float(v[mask].mean()) if mask.any() else None for k,v in errors.items()}}
 result={'n_molecules':len(out),'missing_label_rows':missing,'metrics':metrics,'paired':paired,'subsets':subsets,'scope':'External transfer evaluation; solvent-volume conversion approximation, pH/solid-form uncertainty and source dependence limit comparability. Fixed five-seed average, not TDC five-run metric.'};write(r/'external_results.json',result);write(r/'external_complete.json',{'results_sha256':h(r/'external_results.json')});print(json.dumps(result),flush=True)
if __name__=='__main__':main()
