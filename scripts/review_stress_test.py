"""Display saved stress-test proposals and final matches; no inference or refit."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np
from review_two_slice_validation import display_indices,spread_annotations


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=ROOT/'results/validation/rotation_mismatch_stress_v1')
    p.add_argument('--input',type=Path)
    args=p.parse_args();study=args.study;root=args.input or study/'matching'
    plan=json.loads((study/'study_plan.json').read_text());q=plan['parameters']
    reports=json.loads((root/'reports.json').read_text())
    for report in reports:
        case=report['case'];control=case=='histology_control';out=root/case
        with np.load(out/'predictions.npz') as d:pred={k:d[k] for k in d.files}
        with np.load(study/'pair/inputs.npz') as d:images={k:d[k] for k in d.files}
        if control:
            with np.load(ROOT/plan['source_pair']/'images.npz') as d:images['mri']=d['histology'];images['mri_valid']=d['histology_valid']
        fig,axs=plt.subplots(2,2,figsize=(14,11))
        fig.subplots_adjust(left=.045,right=.98,top=.86,bottom=.17,hspace=.34,wspace=.20)
        for row,(key,errors,metrics,label) in enumerate([
                ('input_pairs',report['input_errors_mm'],report['input'],'INPUT top-3 proposals'),
                ('retained_pairs',report['retained_errors_mm'],report['retained'],'FINAL one-to-one matches')]):
            pairs=pred[key];display=display_indices(pairs,pred['similarity'])
            for ax,name,title in zip(axs[row],['mri','histology'],['Original histology' if control else 'MRI','Transformed histology']):
                im=images[name];valid=images[name+'_valid'];lo,hi=np.percentile(im[valid],[1,99])
                ax.imshow(np.ma.array(im,mask=~valid),cmap='gray',vmin=lo,vmax=hi)
                ax.set_facecolor('#ddd');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f'{title}: {label}',fontsize=11)
            for k in display:
                i,j=pairs[k];error=errors[k];color='#168a35' if error<=2 else '#d28b00' if error<=5 else '#d12c3f'
                for ax,point in zip(axs[row],[pred['source_candidates'][i],pred['target_candidates'][j]]):
                    ax.plot(*point,'o',color=color,ms=4)
                    ax.annotate(str(k+1),point,xytext=(4,4),textcoords='offset points',fontsize=8,color=color,
                        bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=.4))
                fig.add_artist(ConnectionPatch(xyA=pred['source_candidates'][i],coordsA=axs[row,0].transData,
                    xyB=pred['target_candidates'][j],coordsB=axs[row,1].transData,color=color,lw=.8,alpha=.55,clip_on=False))
            axs[row,0].text(.5,-.10,f'ALL {metrics["count"]} pairs: {metrics["count_within_2mm"]} ≤2 mm; {metrics["count_within_5mm"]} ≤5 mm.',transform=axs[row,0].transAxes,ha='center',fontsize=10)
            axs[row,1].text(.5,-.10,f'Showing {len(display)} score-selected pairs with distinct endpoints.',transform=axs[row,1].transAxes,ha='center',fontsize=10)
            if not len(pairs):axs[row,1].text(.5,.5,'No model: insufficient consensus',transform=axs[row,1].transAxes,ha='center',bbox=dict(facecolor='white',alpha=.9))
        rec=report['transform_recovery'];err=f'{rec["median_whole_candidate_transform_error_mm"]:.2f} mm' if rec else 'no model'
        label='HISTOLOGY CONTROL' if control else 'MRI–HISTOLOGY STRESS CASE'
        fig.suptitle(f'{label}: {q["angle"]:+g}° CW, {q["scale"]:.2f}× scale, shift {q["shift"]} px\nEVALUATION colours: green ≤2 mm; amber 2–5 mm; red >5 mm. Median candidate-field recovery error: {err}',fontsize=12,y=.97)
        fig.text(.5,.03,'Top: many-to-many descriptor proposals. Bottom: algorithm-selected one-to-one matches. Labels index pairs within each row.\nUp to 12 highest-score pairs with distinct endpoints shown for readability; all pairs contribute to metrics. Cosine is not confidence.\nFrozen descriptor: radii 16/24/36 px, 41×41 samples, local SSD 5×5; rotation search every 30°. No planted parameters enter inference.\nReference transforms used ONLY after BOTH predictions saved. Approximate alignment; one stress case, not worst-case proof.\nBigBrain — Amunts et al., Science (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=9)
        spread_annotations(fig);fig.savefig(out/'matches.png',dpi=120);plt.close(fig)
    print('Saved stress/matches.png and histology_control/matches.png')

if __name__=='__main__':main()
