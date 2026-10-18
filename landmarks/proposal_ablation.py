"""Fixed proposal-filter ablation on image-derived similarity matrices only."""
import numpy as np
from landmarks.matching import propose

VARIANTS = ('current', 'similarity_only', 'all_mutual')
CONFIG = dict(variants=dict(current='Mutual + cosine >=0.75 + distance ratio <0.9',
    similarity_only='Mutual + cosine >=0.75; no ratio test',
    all_mutual='All valid mutual nearest neighbours; no score or ratio cutoff'),
    ransac=dict(residual_target_px=6.,minimum_consensus=4,max_trials=5000,seed=0),
    evaluation_thresholds_mm=[2.,5.],
    scope='Fixed ablation of saved image-only proposals; no descriptor/detector/RANSAC changes')


def select_variants(scores):
    """Reuse original selection exactly; remove conditions only in named variants."""
    scores=np.asarray(scores)
    mutual,current,_=propose(scores,max_ratio=.9,min_similarity=.75)
    selected=dict(current=current,
        similarity_only=mutual[scores[mutual[:,0],mutual[:,1]]>=.75],
        all_mutual=mutual.copy())
    distance=np.where(scores>=0,np.sqrt(np.maximum(0,2-2*scores)),np.inf)
    audit=[]
    for i,j in mutual:
        alternatives=np.delete(distance[i],j)
        second=float(np.min(alternatives)) if len(alternatives) else np.inf
        defined=bool(np.isfinite(second) and second>1e-12)
        ratio=float(distance[i,j]/second) if defined else None
        audit.append(dict(source_index=int(i),target_index=int(j),similarity=float(scores[i,j]),
            distance_ratio=ratio,ratio_defined=defined,
            passes_ratio=defined and ratio<.9,passes_similarity=bool(scores[i,j]>=.75)))
    return mutual,selected,audit
