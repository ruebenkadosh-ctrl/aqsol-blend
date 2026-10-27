from pathlib import Path
import json,hashlib,requests,pandas as pd
from rdkit import Chem,RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize as st
RDLogger.DisableLog('rdApp.*');r=Path(__file__).resolve().parent
meta=json.loads((r/'source_schema.json').read_text());raw=r/'BigSolDBv2.2.csv'
if not raw.exists():
 response=requests.get(meta['url'],timeout=60);response.raise_for_status();raw.write_bytes(response.content)
f=pd.read_csv(raw,usecols=['SMILES_Solute','Temperature_K','Solvent','SMILES_Solvent','Source','solid_form_label'])
# Eligibility uses metadata/structures only. Outcomes are not loaded.
water=f[f.SMILES_Solvent.eq('O')].copy();near=water[water.Temperature_K.between(297.15,299.15)].copy()
failures=[]
def identity(s):
 m=Chem.MolFromSmiles(s)
 if m is None:return None
 try:
  m=st.Cleanup(m);m=st.FragmentParent(m);m=st.Uncharger().uncharge(m);m=st.TautomerEnumerator().Canonicalize(m);Chem.RemoveStereochemistry(m)
 except Exception as e:
  failures.append({"smiles":s,"error":type(e).__name__});return None
 return Chem.MolToSmiles(m)
a=pd.concat([pd.read_csv('work/tdc/data/admet_group/solubility_aqsoldb/'+x+'.csv',usecols=['Drug']) for x in ['train_val','test']])
lookup={s:identity(s) for s in set(a.Drug)|set(near.SMILES_Solute)};known={lookup[s] for s in a.Drug};near['parent_identity']=near.SMILES_Solute.map(lookup);near['eligible']=near.parent_identity.notna()&~near.parent_identity.isin(known)
near.to_csv(r/'structure_eligibility.csv',index_label='source_row')
result={'dataset_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'total_rows':len(f),'water_rows':len(water),'water_297_15_to_299_15K_rows':len(near),'eligible_rows':int(near.eligible.sum()),'eligible_unique_parent_identities':int(near.loc[near.eligible,'parent_identity'].nunique()),'aqsoldb_rows_excluded_against':len(a),'standardization_failures':len(failures),'unresolved_aqsoldb_structures':sum(lookup[s] is None for s in a.Drug),'endpoint_warning':'Aqueous logS may differ by pH, solid form and conversion; compatibility review pending. No error metrics calculated.'}
(r/'feasibility.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

(r/'standardization_failures.json').write_text(json.dumps(failures,indent=2))
