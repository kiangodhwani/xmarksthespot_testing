"""Descriptor matching and geometric rejection, without access to reference truth."""
import numpy as np
from skimage.measure import ransac
from skimage.transform import SimilarityTransform


def propose(scores,max_ratio=.9,min_similarity=.75):
    """Mutual nearest neighbours, then source-side Euclidean distance ratio.

    Max-over-bank cosine scores are converted to sqrt(2-2*cosine) distances.
    This ratio is a baseline heuristic, not a calibrated match confidence.
    """
    scores=np.asarray(scores)
    if scores.ndim!=2: raise ValueError('Expected pairwise score matrix')
    if not 0<max_ratio<1: raise ValueError('Ratio must be between zero and one')
    if not all(scores.shape): return np.empty((0,2),int),np.empty((0,2),int),np.empty(0)
    distance=np.where(scores>=0,np.sqrt(np.maximum(0,2-2*scores)),np.inf)
    nearest=np.argmin(distance,axis=1);reverse=np.argmin(distance,axis=0)
    mutual=[];accepted=[];ratios=[]
    for i,j in enumerate(nearest):
        if not np.isfinite(distance[i,j]) or reverse[j]!=i: continue
        mutual.append((i,j))
        alternatives=np.delete(distance[i],j)
        second=np.min(alternatives) if len(alternatives) else np.inf
        # At least two usable alternatives are necessary to establish distinctiveness.
        if not np.isfinite(second) or second<=1e-12: continue
        ratio=distance[i,j]/second
        if ratio<max_ratio and scores[i,j]>=min_similarity:
            accepted.append((i,j));ratios.append(float(ratio))
    return np.asarray(mutual,int).reshape(-1,2),np.asarray(accepted,int).reshape(-1,2),np.asarray(ratios)


def noncollinear(a,b):
    for p in [a,b]:
        singular=np.linalg.svd(p-p.mean(axis=0),compute_uv=False)
        if len(singular)<2 or singular[0]<1e-8 or singular[1]/singular[0]<.02: return False
    return True


def geometric_filter(source,target,residual_px=6.,seed=0):
    """Estimate MRI -> histology similarity. Require >=4 noncollinear supporters."""
    source=np.asarray(source,float);target=np.asarray(target,float)
    keep=np.zeros(len(source),bool)
    diagnostic=dict(status='not_enough_proposals',minimum_consensus=4,residual_threshold_histology_px=residual_px,
        seed=seed,max_trials=5000,raw_consensus=0)
    if len(source)<4: return None,keep,diagnostic
    if not noncollinear(source,target):
        diagnostic['status']='degenerate_configuration';return None,keep,diagnostic
    model,inliers=ransac((source,target),SimilarityTransform,min_samples=3,
        residual_threshold=residual_px,is_data_valid=noncollinear,max_trials=5000,
        stop_probability=.999,rng=np.random.default_rng(seed))
    count=int(inliers.sum()) if inliers is not None else 0
    diagnostic['raw_consensus']=count
    if model is None or count<4 or not noncollinear(source[inliers],target[inliers]):
        diagnostic['status']='insufficient_consensus';return None,keep,diagnostic
    # Check residuals against the returned refitted model as well.
    keep=model.residuals(source,target)<residual_px
    if keep.sum()<4:
        diagnostic['status']='insufficient_refit_consensus';return None,np.zeros(len(source),bool),diagnostic
    diagnostic.update(status='estimated_model_not_anatomically_verified',retained=int(keep.sum()),
        scale=float(model.scale),clockwise_degrees=float(np.rad2deg(model.rotation)),
        translation_col_row_px=model.translation.tolist(),matrix=model.params.tolist())
    return model.params,keep,diagnostic
