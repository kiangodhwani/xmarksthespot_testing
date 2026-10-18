"""Explain existing descriptor test examples with full-image context.

No new detection or matching. Boxes and crops use the actual radius combination
that produced each displayed score. Reference positions remain evaluation-only.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,ConnectionPatch
import numpy as np
from scipy.ndimage import map_coordinates
from landmarks.descriptors import describe,RADII

OUT=ROOT/'results/descriptors/polar_gradient_y8'
COLORS=['#0077cc','#ce5500']


def best_settings(a,av,b,bv):
    best=(-2.,0,0,0)
    for i in range(len(RADII)):
        for j in range(len(RADII)):
            if not (av[i] and bv[j]): continue
            for shift in range(12):
                score=float(np.sum(a[i]*np.roll(b[j],shift,axis=1)))
                if score>best[0]: best=(score,i,j,shift)
    return best


def crop(image,point,radius):
    rr,cc=np.indices((81,81),dtype=float)
    coords=np.stack([point[1]+(rr-40)*radius/40,point[0]+(cc-40)*radius/40])
    return map_coordinates(image,coords,order=1,mode='nearest')


def main():
    with np.load(ROOT/'results/transforms/coronal_y8_demo/inputs.npz') as d:
        mri=d['mri'];mvalid=d['mri_valid'];transformed=d['histology'];tvalid=d['histology_valid']
    with np.load(ROOT/'results/transforms/coronal_y8_demo/truth/original_histology.npz') as d:
        hist=d['image'];hvalid=d['valid']
    with np.load(OUT/'mri_descriptors.npz') as d:
        mpoints=d['points_col_row'];mbank=d['descriptors'];musable=d['usable']
    summaries=[]
    for target,himage,hmask in [('histology_original',hist,hvalid),('histology_transformed',transformed,tvalid)]:
        key='mri__'+target
        with np.load(OUT/(key+'_evaluation.npz')) as d:
            indices=d['source_candidate_indices'];positions=d['evaluation_target_points']
            scores=d['similarity'];negatives=d['negative_eligibility']
        eligible=np.flatnonzero(negatives.sum(axis=1)>0)
        query=int(eligible[len(eligible)//2]);negative=int(np.argmax(np.where(negatives[query],scores[query],-2)))
        candidate=int(indices[query]);source=mpoints[candidate]
        partners=positions[[query,negative]]
        hbank,hu=describe(himage,hmask,partners)
        settings=[best_settings(mbank[candidate],musable[candidate],hbank[k],hu[k]) for k in range(2)]
        for k,j in enumerate([query,negative]):
            if not np.isclose(settings[k][0],scores[query,j],atol=1e-6): raise ValueError('Displayed score does not match saved evaluation')
        fig,axs=plt.subplots(3,2,figsize=(12,16),gridspec_kw={'height_ratios':[1.45,1,1]})
        for ax,image,valid,title in [(axs[0,0],mri,mvalid,'1. Full MRI — one candidate centre'),
                                     (axs[0,1],himage,hmask,'2. Full histology — two locations')]:
            lo,hi=np.percentile(image[valid],[1,99])
            ax.imshow(np.ma.array(image,mask=~valid),cmap='gray',vmin=lo,vmax=hi,interpolation='nearest')
            ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.set_facecolor('#ddd')
        for k,(name,label) in enumerate([('reference','B'),('different_location','C')]):
            score,si,ti,shift=settings[k];sr=RADII[si];tr=RADII[ti];color=COLORS[k]
            for ax,point,radius,text in [(axs[0,0],source,sr,f'A{k+1}'),(axs[0,1],partners[k],tr,label)]:
                ax.add_patch(Rectangle((point[0]-radius,point[1]-radius),2*radius,2*radius,fill=False,edgecolor=color,lw=2))
                ax.plot(*point,'+',color=color,ms=9,mew=1.5)
                dx=25
                dy=-45 if k==0 else 45
                ax.annotate(text,point,xytext=(dx,dy),textcoords='offset points',fontsize=12,color=color,
                    bbox=dict(facecolor='white',edgecolor=color,alpha=.95),arrowprops=dict(arrowstyle='->',color=color))
            for ax,image,point,radius,title in [
                (axs[k+1,0],mri,source,sr,f'A{k+1}: MRI crop — radius {sr:g} px'),
                (axs[k+1,1],himage,partners[k],tr,f'{label}: '+('Expected location from alignment' if k==0 else 'Different location preferred by descriptor')+f'\nSimilarity {score:.3f} — radius {tr:g} px')]:
                # Modality-wide display limits keep contrast consistent with full image.
                valid=mvalid if image is mri else hmask
                lo,hi=np.percentile(image[valid],[1,99])
                ax.imshow(crop(image,point,radius),cmap='gray',vmin=lo,vmax=hi,interpolation='nearest')
                ax.axvline(40,color=color,alpha=.5,lw=.6);ax.axhline(40,color=color,alpha=.5,lw=.6)
                ax.set_title(title,color=color,fontsize=11);ax.set_xticks([]);ax.set_yticks([])
                for spine in ax.spines.values(): spine.set_color(color);spine.set_linewidth(2)
        success=settings[0][0]>settings[1][0]
        outcome='Reference patch scores higher.' if success else 'Failure: the different location scores higher.'
        fig.suptitle('Read top to bottom: full images → boxed crops → descriptor preference\n'+outcome,fontsize=15)
        fig.subplots_adjust(left=.04,right=.97,top=.90,bottom=.15,hspace=.32,wspace=.15)
        fig.text(.5,.054,'A1 and A2 share the SAME MRI centre; patch sizes may differ because the comparison searches scales.\nB is expected from the approximate alignment. C is >5 mm away; neither is an accepted predicted match.\nCrops retain native orientation. Scores search all radii and 30° angular shifts; they are not probabilities.',ha='center',fontsize=9)
        fig.text(.5,.014,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=9)
        output=OUT/('context_'+target+'.png')
        fig.savefig(output,dpi=150);plt.close(fig)
        summaries.append(dict(file=output.name,source_candidate=candidate,source_col_row=source.tolist(),
            reference_col_row=partners[0].tolist(),different_location_col_row=partners[1].tolist(),
            settings=[dict(score=v[0],source_radius_px=RADII[v[1]],target_radius_px=RADII[v[2]],sector_roll=v[3]) for v in settings]))
        print(output)
    (OUT/'context_examples.json').write_text(json.dumps(summaries,indent=2))

if __name__=='__main__': main()
