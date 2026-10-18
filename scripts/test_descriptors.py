"""Compute candidate signatures independently; evaluate known-location patches.

Descriptor validation only. Evaluation locations are supplied by saved geometry,
NOT predicted matches. Spatially ambiguous negatives within 5 mm are excluded.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import hashlib
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial.distance import cdist
from landmarks.descriptors import describe,similarity,RADII
from landmarks.transforms import apply_points


def example_patch(image,point,radius=24):
    y,x=np.indices((65,65),dtype=float); y=(y-32)*radius/32+point[1];x=(x-32)*radius/32+point[0]
    return map_coordinates(image,np.stack([y,x]),order=1,mode='nearest')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/descriptors/polar_gradient_y8')
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Choose a new --output directory')
    input_path=ROOT/'results/transforms/coronal_y8_demo/inputs.npz'
    digest=hashlib.sha256(input_path.read_bytes()).hexdigest()
    with np.load(input_path) as d:
        images={'mri':d['mri'],'histology_transformed':d['histology']}
        valid={'mri':d['mri_valid'],'histology_transformed':d['histology_valid']}
    with np.load(ROOT/'results/transforms/coronal_y8_demo/truth/original_histology.npz') as d:
        images['histology_original']=d['image'];valid['histology_original']=d['valid']
    points={name:np.loadtxt(ROOT/f'results/candidates/harris_y8_peakfix/{name}.csv',delimiter=',',skiprows=1)[:,:2] for name in images}
    args.output.mkdir(parents=True)
    descriptor_stats={};descriptors={}
    for name in images:
        bank,usable=describe(images[name],valid[name],points[name]);descriptors[name]=(bank,usable)
        np.savez_compressed(args.output/(name+'_descriptors.npz'),points_col_row=points[name],descriptors=bank,usable=usable,radii_px=RADII)
        descriptor_stats[name]=dict(candidates=len(bank),with_usable_signature=int(usable.any(axis=1).sum()),usable_per_radius=usable.sum(axis=0).tolist())
    # Only AFTER independent signatures are saved, read evaluation geometry.
    geometry=json.loads((ROOT/'results/transforms/coronal_y8_demo/truth/geometry.json').read_text())
    t=np.asarray(geometry['original_to_transformed'])
    cases=[('histology_original','histology_transformed',t),('mri','histology_original',np.eye(3)),('mri','histology_transformed',t)]
    evaluation={};distributions=[];examples=[]
    for a,b,matrix in cases:
        projected=apply_points(points[a],matrix)
        target,tvalid=describe(images[b],valid[b],projected)
        source,svalid=descriptors[a]
        keep=svalid.any(axis=1)&tvalid.any(axis=1)
        ids=np.flatnonzero(keep)
        src=source[keep];sv=svalid[keep];target=target[keep];tv=tvalid[keep]
        scores=similarity(src,sv,target,tv)
        source_points=points[a][keep];target_points=projected[keep]
        distances=cdist(source_points*.2,source_points*.2)
        negatives=distances>5.
        diagonal=np.diag(scores)
        eligible=negatives.sum(axis=1)>0
        ranks=1+((scores>=diagonal[:,None]-1e-7)&negatives).sum(axis=1) # ties conservatively count against retrieval
        truewins=[]
        for i in np.flatnonzero(eligible):
            neg=scores[i,negatives[i]]
            truewins.append(float(np.mean((diagonal[i]>neg+1e-7)+.5*(np.abs(diagonal[i]-neg)<=1e-7))))
        name=a+'__'+b
        evaluation[name]=dict(queries=int(eligible.sum()),excluded_unusable=int((~keep).sum()),
            true_partner_top1_fraction=float(np.mean(ranks[eligible]==1)),
            mean_random_top1_chance=float(np.mean(1/(negatives.sum(axis=1)[eligible]+1))),
            median_true_partner_rank=float(np.median(ranks[eligible])),
            mean_true_vs_negative_win_fraction=float(np.mean(truewins)),
            median_positive_similarity=float(np.median(diagonal[eligible])),
            median_negative_similarity=float(np.median(scores[negatives])),
            median_number_of_negatives=float(np.median(negatives.sum(axis=1)[eligible])),
            caveat='Controlled known-location gallery, not candidate matching: target patch locations supplied by reference geometry; negatives within 5 mm excluded')
        distributions.append((name,diagonal[eligible],scores[negatives]))
        # A deterministic illustrative query: median source CSV index, not best result.
        i=int(np.flatnonzero(eligible)[len(np.flatnonzero(eligible))//2])
        j=int(np.argmax(np.where(negatives[i],scores[i],-2)))
        examples.append((name,example_patch(images[a],source_points[i]),example_patch(images[b],target_points[i]),
            example_patch(images[b],target_points[j]),float(diagonal[i]),float(scores[i,j]),int(ids[i])))
        np.savez_compressed(args.output/(name+'_evaluation.npz'),source_candidate_indices=ids,
            evaluation_target_points=target_points,similarity=scores,negative_eligibility=negatives)
        print(name,evaluation[name],flush=True)
    report=dict(method='Unsigned polar-gradient histogram, square patch resampling then four rings x twelve sectors x eight relative gradient directions',
        descriptor_dimensions_per_radius=384,radii_pixels=list(RADII),
        comparison='Maximum cosine over 3x3 radius combinations and 12 discrete angular shifts; parameters fixed, not selected from planted transform',
        rotation_resolution_degrees=30,scale_support='Discrete patch-size bank, not guaranteed continuous scale invariance',
        input_sha256=digest,descriptor_stats=descriptor_stats,evaluation=evaluation,
        limitations=['Known-location test isolates descriptor behaviour from detector localisation',
            'MRI/histology positives are approximate reference correspondences',
            'Nearby negatives within 5 mm excluded, so local ambiguity remains untested',
            'One slice and one transform; no held-out validation',
            'Edge-direction invariance to global inversion does not make this modality independent'],
        attribution=geometry['attribution'])
    (args.output/'report.json').write_text(json.dumps(report,indent=2))
    fig,axs=plt.subplots(3,3,figsize=(11,10),layout='constrained')
    for row,(name,source,positive,negative,ps,ns,index) in enumerate(examples):
        for col,(patch,title) in enumerate([(source,f'Source candidate #{index}'),(positive,f'Reference partner\nsimilarity {ps:.3f}'),(negative,f'Most similar negative (>5 mm away)\nsimilarity {ns:.3f}')]):
            axs[row,col].imshow(patch,cmap='gray',interpolation='nearest');axs[row,col].set_title(title,fontsize=10)
            axs[row,col].axvline(32,color='#e60070',alpha=.25,lw=.5);axs[row,col].axhline(32,color='#e60070',alpha=.25,lw=.5)
            axs[row,col].set_xticks([]);axs[row,col].set_yticks([])
        axs[row,0].set_ylabel(name.replace('__',' → ').replace('_',' '),fontsize=8)
    fig.suptitle('What the descriptor compares — 24 px-radius examples\nReference locations are evaluation inputs, not predicted matches; comparison uses all patch sizes')
    fig.get_layout_engine().set(rect=(0,.05,1,.95));fig.text(.5,.015,geometry['attribution'],ha='center',fontsize=8)
    fig.savefig(args.output/'patch_examples.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    for ax,(name,pos,neg) in zip(axs,distributions):
        ax.hist(neg,bins=30,density=True,alpha=.5,label='Unrelated patches (>5mm)')
        ax.hist(pos,bins=20,density=True,alpha=.5,label='Reference partners')
        ax.set_title(name.replace('__',' → ').replace('_',' '),fontsize=10)
        ax.set_xlabel('Signature similarity (higher = more alike)');ax.set_ylabel('Density');ax.legend(fontsize=7)
    fig.suptitle('Descriptor-only validation; distributions can overlap')
    fig.savefig(args.output/'similarities.png',dpi=140);plt.close(fig)
    assert hashlib.sha256(input_path.read_bytes()).hexdigest()==digest
    print('Saved:',args.output)

if __name__=='__main__': main()
