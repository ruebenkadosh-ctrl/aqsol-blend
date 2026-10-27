from pathlib import Path
from collections import Counter
import json,pandas as pd,numpy as np
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*');r=Path(__file__).resolve().parent
f=pd.read_csv(r/'structure_eligibility.csv');fail=json.loads((r/'standardization_failures.json').read_text())
def heavy_signature(s):
 m=Chem.MolFromSmiles(s)
 if m is None:return None
 fragments=Chem.GetMolFrags(m,asMols=True);m=max(fragments,key=lambda x:x.GetNumHeavyAtoms())
 return tuple(sorted(Counter(a.GetAtomicNum() for a in m.GetAtoms() if a.GetAtomicNum()>1).items()))
sigs={heavy_signature(x['smiles']) for x in fail};f['fallback_excluded']=f.SMILES_Solute.map(lambda s:heavy_signature(s) in sigs)
f['eligible_final']=f.eligible&~f.fallback_excluded
# Over-exclude every matching parent heavy-element composition, regardless of charge/tautomer/bonds.
e=f[f.eligible_final].copy();dev=pd.read_csv('work/tdc/data/admet_group/solubility_aqsoldb/train_val.csv',usecols=['Drug'])
def fp(s):return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s),2,nBits=2048)
trainfp=[fp(s) for s in dev.Drug];sim={s:max(DataStructs.BulkTanimotoSimilarity(fp(s),trainfp)) for s in e.SMILES_Solute.unique()};e['max_train_tanimoto']=e.SMILES_Solute.map(sim)
f.to_csv(r/'exclusion_ledger.csv',index=False);e.to_csv(r/'eligible_metadata.csv',index=False)
q=e.drop_duplicates('parent_identity');out={'rows':len(e),'unique_parent_identities':len(q),'extra_exclusions_due_to_failed_reference_structure':int((f.eligible&f.fallback_excluded).sum()),'fallback':'Conservative exclusion by largest-fragment heavy-element composition; independent of bond order, hydrogen, charge and stereochemistry. Not exact identity resolution.','nearest_training_similarity_quantiles':np.quantile(q.max_train_tanimoto,[0,.25,.5,.75,1]).tolist(),'temperature_minmax':[float(e.Temperature_K.min()),float(e.Temperature_K.max())],'missing_solid_form_rows':int(e.solid_form_label.isna().sum()),'pH_metadata':'not provided in source schema','labels_loaded':False}
(r/'cohort_review.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
