from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
s=json.loads((ROOT/'analysis/summary.json').read_text(encoding='utf-8'))
keys=['RelativeError','EventAmplitudeError','WeakEventError','AVOGradientError','NoiseWindowRMS','FaultNeighborhoodError','PinchoutTailError']
labels=['Waveform error','Event amplitude error','Weak event error','AVO gradient error','Noise window RMS','Fault neighborhood error','Pinchout tail error']
fig,ax=plt.subplots(figsize=(9,5.6));y=np.arange(len(keys))
profiles=[('SourceMixture',-.16,'#0072B2','o'),('SourceOnly',0,'#D55E00','s')]
if 'Consensus' in s['means']:profiles.append(('Consensus',.16,'#009E73','^'))
for method,dy,color,marker in profiles:
    ratio=[s['means'][method][k]/s['means']['Baseline'][k] for k in keys]
    ax.scatter(ratio,y+dy,c=color,marker=marker,s=45,label=method,zorder=3)
ax.axvline(1,color='#555555',ls='--',lw=1,label='Current baseline')
ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(left=0)
ax.set_xlabel('Ratio of mean metric to baseline (lower is better)')
ax.set_title('One fixed rule across all 60 paired development cases',loc='left',fontsize=12)
ax.grid(axis='x',alpha=.2);ax.spines[['top','right']].set_visible(False)
ax.legend(loc='best',fontsize=9)
fig.text(.01,.01,'Each metric retains its original scoring subset; weak: 12 cases, AVO: 24 cases. No per-model selection.',fontsize=8)
fig.tight_layout(rect=[0,.04,1,1]);out=ROOT/'figures';out.mkdir(exist_ok=True)
fig.savefig(out/'global_candidate_metrics.png',dpi=200);fig.savefig(out/'global_candidate_metrics.pdf')
