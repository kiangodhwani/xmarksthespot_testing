"""Render frozen new-slice predictions and separate reference evaluation."""
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
from matplotlib.patches import ConnectionPatch,Patch
import numpy as np


def display_indices(pairs,scores,limit=12):
    """Score-only display selection; distinct endpoints to avoid stacked labels."""
    if not len(pairs):return []
    order=np.argsort(-scores[pairs[:,0],pairs[:,1]],kind='stable')
    selected=[];sources=set();targets=set()
    for index in order:
        i,j=pairs[index]
        if i in sources or j in targets:continue
        selected.append(int(index));sources.add(i);targets.add(j)
        if len(selected)==limit:break
    return selected


def spread_annotations(fig):
    """Keep numeric labels readable when detected points are close."""
    from matplotlib.text import Annotation
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    offsets=[(4,4),(4,-11),(-12,4),(-12,-11),(10,11),(-18,11),(10,-18),(-18,-18),(19,0),(-24,0)]
    for ax in fig.axes:
        placed=[]
        for label in [t for t in ax.texts if isinstance(t,Annotation)]:
            for offset in offsets:
                label.set_position(offset);bounds=label.get_window_extent(renderer).expanded(1.12,1.18)
                if not any(bounds.overlaps(other) for other in placed):break
            placed.append(bounds)


def case_figure(study,root,slice_id,case,report):
    folder=study/slice_id;control=case=='histology_control'
    pair_dir=folder/'fixed' if control or case=='fixed' else folder/'random_seed42'/case
    out=root/slice_id/case
    with np.load(out/'predictions.npz') as d:pred={k:d[k] for k in d.files}
    with np.load(pair_dir/'inputs.npz') as d:images={k:d[k] for k in d.files}
    if control:
        with np.load(folder/'pair/images.npz') as d:images['mri']=d['histology'];images['mri_valid']=d['histology_valid']
    fig,axs=plt.subplots(2,2,figsize=(14,11))
    fig.subplots_adjust(left=.045,right=.98,top=.87,bottom=.16,hspace=.36,wspace=.21)
    for row,(key,errors,stage,title) in enumerate([
            ('input_pairs',report['input_errors_mm'],report['input'],'INPUT top-3 proposals'),
            ('retained_pairs',report['retained_errors_mm'],report['retained'],'FINAL one-to-one matches')]):
        pairs=pred[key];show=display_indices(pairs,pred['similarity'])
        for ax,name,label in zip(axs[row],['mri','histology'],['Original histology' if control else 'MRI','Transformed histology']):
            image=images[name];valid=images[name+'_valid'];lo,hi=np.percentile(image[valid],[1,99])
            ax.imshow(np.ma.array(image,mask=~valid),cmap='gray',vmin=lo,vmax=hi)
            ax.set_facecolor('#ddd');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f'{label}: {title}',fontsize=11)
        for k in show:
            i,j=pairs[k];error=errors[k]
            color='#168a35' if error<=2 else '#d28b00' if error<=5 else '#d12c3f'
            for ax,point in zip(axs[row],[pred['source_candidates'][i],pred['target_candidates'][j]]):
                ax.plot(*point,'o',color=color,ms=4)
                ax.annotate(str(k+1),point,xytext=(4,4),textcoords='offset points',fontsize=8,color=color,
                    bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=.4))
            fig.add_artist(ConnectionPatch(xyA=pred['source_candidates'][i],coordsA=axs[row,0].transData,
                xyB=pred['target_candidates'][j],coordsB=axs[row,1].transData,color=color,lw=.8,alpha=.55,clip_on=False))
        axs[row,0].text(.5,-.10,f'ALL {stage["count"]} pairs: {stage["count_within_2mm"]} ≤2 mm; {stage["count_within_5mm"]} ≤5 mm.',transform=axs[row,0].transAxes,ha='center',fontsize=10)
        axs[row,1].text(.5,-.10,f'Showing {len(show)} score-selected pairs with distinct endpoints.',transform=axs[row,1].transAxes,ha='center',fontsize=10)
        if not len(pairs):axs[row,1].text(.5,.5,'No model: insufficient consensus',transform=axs[row,1].transAxes,ha='center',bbox=dict(facecolor='white',alpha=.95))
    y=-20 if slice_id=='y_m20' else 35
    rec=report['transform_recovery'];err=f'{rec["median_whole_candidate_transform_error_mm"]:.2f} mm' if rec else 'no model'
    fig.suptitle(f'New slice y={y} mm — {case}: frozen MIND top-3 + relational geometry\nEVALUATION-only colours: green ≤2 mm; amber 2–5 mm; red >5 mm. Median candidate-field recovery error: {err}',fontsize=12,y=.97)
    fig.text(.5,.028,'Top row is many-to-many proposal input; bottom row is algorithm-selected one-to-one output. Numbers index pairs within each row.\nDisplay uses cosine score only, greedily skipping repeated endpoints for readability; ALL pairs contribute to metrics.\nScores are not probabilities. Descriptor settings unchanged: 3×3 radii, 12 rotation samples; support sizes documented in MIND_EXPERIMENT.md.\nReference transforms used only after ALL new-slice predictions saved. Alignment approximate; same specimen as development.\nBigBrain — Amunts et al., Science (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    spread_annotations(fig)
    fig.savefig(out/'matches.png',dpi=120);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=ROOT/'results/validation/two_slices_v1')
    p.add_argument('--input',type=Path)
    args=p.parse_args();root=args.input or args.study/'matching'
    reports=json.loads((root/'reports.json').read_text())
    for r in reports:case_figure(args.study,root,r['slice_id'],r['case'],r)
    lookup={(r['slice_id'],r['case']):r for r in reports};cases=['fixed']+[f'pair_{i:03d}' for i in range(1,6)]
    fig,axs=plt.subplots(1,2,figsize=(14,6),sharey=True)
    fig.subplots_adjust(left=.06,right=.99,top=.78,bottom=.27,wspace=.16)
    for ax,slice_id,y in zip(axs,['y_m20','y35'],[-20,35]):
        data=[lookup[(slice_id,c)]['retained'] for c in cases]
        green=np.array([r['count_within_2mm'] for r in data]);amber=np.array([r['count_within_5mm'] for r in data])-green
        total=np.array([r['count'] for r in data]);red=total-green-amber;x=np.arange(6)
        ax.bar(x,green,color='#168a35');ax.bar(x,amber,bottom=green,color='#d28b00');ax.bar(x,red,bottom=green+amber,color='#d12c3f')
        for i in x:ax.text(i,total[i]+.2,f'{green[i]}/{total[i]}',ha='center',fontsize=10)
        ax.set_xticks(x,['Fixed','001','002','003','004','005']);ax.set_title(f'New coronal slice y={y} mm');ax.set_ylabel('Final one-to-one matches')
        ax.set_ylim(0,max(r['retained']['count'] for r in reports if r['case']!='histology_control')+3)
    fig.suptitle('Frozen-pipeline validation on two additional slices\nLabels = within 2 mm / total retained; colours use reference only after inference',fontsize=14,y=.98)
    fig.legend(handles=[Patch(color='#168a35',label='≤2 mm'),Patch(color='#d28b00',label='2–5 mm'),Patch(color='#d12c3f',label='>5 mm')],loc='upper center',bbox_to_anchor=(.5,.89),ncol=3)
    control=[]
    for sid,y in [('y_m20',-20),('y35',35)]:
        r=lookup[(sid,'histology_control')]['retained'];control.append(f'y={y}: {r["count_within_2mm"]}/{r["count"]}')
    fig.text(.5,.04,'Histology controls within 2 mm / retained: '+ '; '.join(control)+'.\nSame parameters and six transformations per slice; no tuning. Two additional slices of ONE specimen, not new subjects.\nReturned model ≠ correct registration: inspect retained errors, spatial coverage and whole-candidate recovery.\nBigBrain — Amunts et al., Science (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=10)
    fig.savefig(root/'comparison.png',dpi=140);plt.close(fig)
    groups=[]
    for sid in ['y_m20','y35']:
        rs=[r for r in reports if r['slice_id']==sid and r['case']!='histology_control']
        groups.append(dict(slice_id=sid,cases=len(rs),models_returned=sum(r['retained']['count']>0 for r in rs),
            retained=sum(r['retained']['count'] for r in rs),within_2mm=sum(r['retained']['count_within_2mm'] for r in rs),
            within_5mm=sum(r['retained']['count_within_5mm'] for r in rs),
            recovery_median_error_mm_by_case={r['case']:r['transform_recovery']['median_whole_candidate_transform_error_mm'] if r['transform_recovery'] else None for r in rs}))
    (root/'slice_summary.json').write_text(json.dumps(groups,indent=2))
    print(json.dumps(groups,indent=2))
    print('Saved 14 full-image match figures and all-case comparison')

if __name__=='__main__':main()
