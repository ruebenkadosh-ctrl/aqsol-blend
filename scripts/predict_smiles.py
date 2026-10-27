"""Predict with a trusted frozen AqSol-Blend package; no label access."""
import argparse,json,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--models',type=Path,required=True);p.add_argument('--smiles',nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if a.output.exists():p.error('Output exists; choose a new file')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from predict import predict
out=predict(a.smiles,model_dir=a.models,cache_dir=a.output.parent/'cache')
a.output.write_text(json.dumps({'smiles':a.smiles,'per_seed':out['per_seed'].tolist(),'mean':out['mean'].tolist()},indent=2))
