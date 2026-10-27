import os
os.environ['OMP_NUM_THREADS']='6'
from pathlib import Path
import sys,json,hashlib
import pandas as pd,numpy as np,joblib
r=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
with (out/'rf_started.json').open('x') as f:json.dump({'scope':'Post hoc single fixed RF test baseline; no ensemble rescoring or tuning','rf_sha256':hashlib.sha256((r/'work/aqsol_external/rf.joblib').read_bytes()).hexdigest()},f)
sys.path.insert(0,str(r/'outputs/solubility_mordred_candidate'));from features import featurize
path=r/'work/tdc/data/admet_group/solubility_aqsoldb/test.csv'
smiles=pd.read_csv(path,usecols=['Drug']).Drug.tolist();x,meta=featurize(smiles,out/'cache');model=joblib.load(r/'work/aqsol_external/rf.joblib');pred=model.predict(x[:,meta['groups']['descriptors'][0]:]);assert np.isfinite(pred).all();np.save(out/'rf_test_predictions.npy',pred)
y=pd.read_csv(path,usecols=['Y']).Y.to_numpy();result={'n':len(y),'mae':float(np.mean(abs(y-pred))),'rmse':float(np.sqrt(np.mean((y-pred)**2))),'scope':'Post hoc single RF fit on all 7985 development rows; not TDC five-seed protocol. Fixed model previously evaluated externally; no tuning.'}
(out/'rf_test_results.json').write_text(json.dumps(result,indent=2));print(result)
