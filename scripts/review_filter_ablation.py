"""Summarise frozen filter-ablation results without refitting or selecting models."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import csv
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from landmarks.proposal_ablation import VARIANTS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=ROOT/'results/matching/filter_ablation_v1')
    args=p.parse_args();root=args.input
    reports=json.loads((root/'reports.json').read_text())
    lookup={(r['method'],r['case'],r['variant']):r for r in reports}
    cases=['coronal_y8_demo']+[f'pair_{i:03d}' for i in range(1,6)]
    titles=['Current: score + ratio','Score only; NO ratio','All mutual; NO score or ratio']
    fig,axs=plt.subplots(2,3,figsize=(16,9),sharey=True)
    fig.subplots_adjust(left=.055,right=.985,top=.83,bottom=.16,hspace=.38,wspace=.15)
    for row,method in enumerate(['mind','polar']):
        for col,variant in enumerate(VARIANTS):
            ax=axs[row,col];data=[lookup[(method,c,variant)]['stages']['retained'] for c in cases]
            green=np.array([d['count_within_2mm'] for d in data]);amber=np.array([d['count_within_5mm'] for d in data])-green
            red=np.array([d['count'] for d in data])-green-amber;x=np.arange(6)
            ax.bar(x,green,color='#168a35');ax.bar(x,amber,bottom=green,color='#d28b00');ax.bar(x,red,bottom=green+amber,color='#d12c3f')
            for i,d in enumerate(data):ax.text(i,d['count']+.2,f'{green[i]}/{d["count"]}',ha='center',fontsize=10)
            ax.set_xticks(x,['Fixed','001','002','003','004','005']);ax.set_ylim(0,13)
            ax.set_title(f'{method.upper()} — {titles[col]}',fontsize=11)
            if col==0:ax.set_ylabel('Final RANSAC pairs')
            ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Proposal-filter ablation: same saved descriptors and candidates; same RANSAC\nBar labels = within 2 mm / total retained; all six cases included',y=.98,fontsize=14)
    fig.legend(handles=[Patch(color='#168a35',label='≤2 mm'),Patch(color='#d28b00',label='2–5 mm'),Patch(color='#d12c3f',label='>5 mm')],loc='upper center',bbox_to_anchor=(.5,.91),ncol=3)
    fig.text(.5,.025,'Reference geometry is EVALUATION ONLY. All retained pairs here are within 5 mm; that does not establish anatomical correctness.\nScore-only and all-mutual inputs are identical here: every mutual proposal already has cosine ≥0.75.\nHistology controls within 2 mm: MIND 65 → 96; polar 117 → 127 after removing ratio.\nOne slice/specimen, six augmentations; fixed seed; no held-out validation. BigBrain — Amunts et al. (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=10)
    fig.savefig(root/'comparison.png',dpi=140);plt.close(fig)
    rows=list(csv.DictReader((root/'mind/pair_004/proposal_audit_evaluation.csv').open()))
    fig,ax=plt.subplots(figsize=(10,7));fig.subplots_adjust(left=.1,right=.97,bottom=.25,top=.85)
    for r in rows:
        if not r['distance_ratio']:continue
        ratio=float(r['distance_ratio']);score=float(r['similarity']);error=float(r['reference_error_mm'])
        color='#168a35' if error<=2 else '#d28b00' if error<=5 else '#d12c3f'
        ax.scatter(ratio,score,color=color,s=38,edgecolors='white',linewidth=.4)
    ax.axvline(.9,color='black',linestyle='--',label='Ratio must be <0.9')
    ax.axhline(.75,color='gray',linestyle=':',label='Cosine must be ≥0.75')
    ax.set_xlim(0,1.03);ax.set_ylim(.73,1.01);ax.set_xlabel('Best / second-best descriptor distance (smaller = more distinctive)')
    ax.set_ylabel('Best cosine similarity (not probability)');ax.legend(loc='lower left')
    fig.suptitle('MIND pair 004: why the current filter rejects mutual proposals\nGreen ≤2 mm; amber 2–5 mm; red >5 mm. All 32 mutual proposals shown.',fontsize=12)
    fig.text(.5,.035,'11 reference-correct mutual pairs: 5 pass; 6 fail the ratio condition. None fail the cosine cutoff.\nEvery point is one mutual proposal, including proposals not shown in the earlier top-12 match view.\nDistances = sqrt(2−2×cosine); colour uses reference only after inference.\nMIND-style configuration is fixed: 3×3 radii and 12 rotation samples; no reference geometry enters scoring.',ha='center',fontsize=10)
    fig.savefig(root/'pair004_filter_audit.png',dpi=140);plt.close(fig)
    print('Saved comparison and pair004 ratio audit')

if __name__=='__main__':main()
