from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
f=pd.read_csv(r/'work/aqsol_external/external_molecule_results.csv');bins=[f.max_train_tanimoto<.3,(f.max_train_tanimoto>=.3)&(f.max_train_tanimoto<.6),f.max_train_tanimoto>=.6];labels=['< 0.3','0.3–< 0.6','≥ 0.6'];result=[]
rng=np.random.default_rng(20261027)
for label,mask in zip(labels,bins):
 sub=f[mask];a=abs(sub.current-sub.truth).to_numpy();b=abs(sub.rf-sub.truth).to_numpy();delta=b-a;idx=rng.integers(0,len(sub),size=(2000,len(sub)))
 result.append({'bin':label,'n':len(sub),'ensemble_mae':float(a.mean()),'rf_mae':float(b.mean()),'gain':float(delta.mean()),'posthoc_paired_bootstrap95':np.quantile(delta[idx].mean(axis=1),[.025,.975]).tolist()})
(out/'similarity_comparison.json').write_text(json.dumps(result,indent=2))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,axes=plt.subplots(1,2,figsize=(8,3.7),sharex=True,sharey=True,layout='constrained')
lo=min(f.truth.min(),f.current.min(),f.rf.min())-.3;hi=max(f.truth.max(),f.current.max(),f.rf.max())+.3
for ax,col,title,color in zip(axes,['current','rf'],['AqSol-Blend','Fixed random forest'],['#176B87','#BD632F']):
 ax.scatter(f.truth,f[col],s=16,alpha=.65,c=color,edgecolors='none');ax.plot([lo,hi],[lo,hi],'--',color='0.45',lw=1);ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Observed log solubility',title=title);ax.set_aspect('equal')
axes[0].set_ylabel('Predicted log solubility');fig.savefig(out/'external_parity.png');plt.close(fig)
fig,ax=plt.subplots(figsize=(6.3,3.8),layout='constrained');x=np.arange(3);w=.34
ax.bar(x-w/2,[v['ensemble_mae'] for v in result],w,label='AqSol-Blend',color='#516A91');ax.bar(x+w/2,[v['rf_mae'] for v in result],w,label='Fixed random forest',color='#B6BFCC');ax.set_xticks(x,[f'{v["bin"]}\n(n={v["n"]})' for v in result]);ax.set(xlabel='Nearest-development Morgan Tanimoto similarity',ylabel='Mean absolute error (log units)',ylim=(0,1.55));ax.legend(frameon=False);fig.savefig(out/'error_similarity.png');plt.close(fig)
print(json.dumps(result,indent=2))
