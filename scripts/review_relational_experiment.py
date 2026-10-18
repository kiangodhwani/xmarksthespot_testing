"""Plot saved proposal/geometry results; never choose or refit a model."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import csv
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=ROOT/'results/matching/relational_v1')
    args=p.parse_args();root=args.input
    with (root/'summary.csv').open() as f:summary=list(csv.DictReader(f))
    with (root/'retrieval.csv').open() as f:retrieval=list(csv.DictReader(f))
    cases=['coronal_y8_demo']+[f'pair_{i:03d}' for i in range(1,6)]
    lookup={(r['case'],r['variant']):r for r in summary};ranks={r['case']:r for r in retrieval}
    fig,axs=plt.subplots(1,3,figsize=(16,6),sharey=True)
    fig.subplots_adjust(left=.05,right=.99,top=.78,bottom=.27,wspace=.12)
    for ax,variant,title in zip(axs,['mutual_ransac','mutual_relational','top3_relational'],
            ['Mutual + seeded RANSAC','Mutual + relational geometry','Top 3 + relational geometry']):
        rows=[lookup[(case,variant)] for case in cases]
        green=np.array([int(r['within_2mm']) for r in rows]);amber=np.array([int(r['within_5mm']) for r in rows])-green
        total=np.array([int(r['retained']) for r in rows]);red=total-green-amber
        x=np.arange(6);ax.bar(x,green,color='#168a35');ax.bar(x,amber,bottom=green,color='#d28b00');ax.bar(x,red,bottom=green+amber,color='#d12c3f')
        for i in x:ax.text(i,total[i]+.2,f'{green[i]}/{total[i]}',ha='center',fontsize=10)
        ax.set_xticks(x,['Fixed','001','002','003','004','005']);ax.set_ylim(0,16);ax.set_title(title,fontsize=11)
    axs[0].set_ylabel('Final one-to-one matches')
    fig.suptitle('MIND: preserve candidate alternatives, then use relational geometry\nBar labels = within 2 mm / total retained; colours are evaluation only',fontsize=14,y=.98)
    fig.legend(handles=[Patch(color='#168a35',label='≤2 mm'),Patch(color='#d28b00',label='2–5 mm'),Patch(color='#d12c3f',label='>5 mm')],loc='upper center',bbox_to_anchor=(.5,.89),ncol=3)
    fig.text(.5,.04,'Descriptors, candidate sets, cosine cutoff 0.75, residual <6 target pixels and minimum four supporters unchanged.\nRelational geometry replaces sampled RANSAC with pair hypotheses and one-to-one assignment; mutual solver control shown.\nTop 3 uses source-side alternatives and removes mutual-nearest restriction. Same slice, six transformations; no held-out validation.\nHistology control within 2 mm: mutual RANSAC 96/96; mutual relational 96/96; top-3 relational 137/137.\nBigBrain — Amunts et al., Science (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(root/'comparison.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(13,6));fig.subplots_adjust(left=.07,right=.98,top=.84,bottom=.27)
    x=np.arange(6)
    for k,offset,color in [(1,-.24,'#0072b2'),(3,0,'#d55e00'),(5,.24,'#6c4591')]:
        values=[int(ranks[c][f'top{k}_sources_with_partner_within_2mm']) for c in cases]
        bars=ax.bar(x+offset,values,.24,label=f'Top {k}',color=color);ax.bar_label(bars,padding=2,fontsize=9)
    available=[int(ranks[c]['sources_with_any_target_within_2mm']) for c in cases]
    ax.plot(x,available,'ko--',label='Any detected target within 2 mm')
    ax.set_ylim(0,50);ax.set_xticks(x,['Fixed','001','002','003','004','005']);ax.set_ylabel('MRI candidates with at least one reference-close target')
    ax.legend(loc='upper right',ncol=2);ax.set_title('Evaluation-only retrieval diagnostic: where do reference-close partners rank?',fontsize=13)
    fig.text(.5,.04,'Each case has 86 MRI candidates. Count a source once if any source-ranked, score-eligible target is within 2 mm.\nTop 1 is source nearest-neighbour retrieval, not mutual matching. Top 5 is diagnostic only; inference uses fixed top 3.\nBlack line is candidate-location availability, not guaranteed anatomical correspondence. No oracle locations inserted into proposals.\nColours distinguish retrieval settings here. Reference alignment approximate; one slice/specimen.',ha='center',fontsize=10)
    fig.savefig(root/'retrieval.png',dpi=140);plt.close(fig)
    print('Saved comparison.png and retrieval.png')

if __name__=='__main__':main()
