"""Image-only alternative proposals and deterministic similarity consensus.

No reference geometry, model scale/range prior or anatomical labels are inputs.
Two distinct correspondences determine a similarity hypothesis. Other proposals
vote via the unchanged <6 target-pixel tolerance; assignments are one-to-one.
"""
from itertools import combinations
import numpy as np
from scipy.optimize import linear_sum_assignment
from skimage.transform import SimilarityTransform
from landmarks.matching import noncollinear

CONFIG=dict(top_k=3,min_similarity=.75,residual_target_px=6.,minimum_consensus=4,
    hypotheses='Exhaustive two-proposal similarity fits; distinct source and target',
    assignment='Maximum cardinality one-to-one among residual<6 proposals, then minimum squared residual',
    selection='Most final assigned pairs, then smallest squared residual, then largest cosine sum',
    refinement='At most three assignment/refit iterations; residual<6; noncollinear on both sides',
    pruning='Rank initial support upper bounds; skip hypotheses below best final count; heuristic after refits, not global optimality guarantee',
    reference_access='None during proposal selection, hypothesis ranking or fit',
    caveat='New deterministic geometry solver; mutual proposals also run through it as a solver control')


def topk_proposals(scores,k=3,min_similarity=.75):
    scores=np.asarray(scores)
    if scores.ndim!=2 or k<1:raise ValueError('Expected score matrix and positive k')
    pairs=[]
    for i,row in enumerate(scores):
        js=np.argsort(-row,kind='stable')[:k]
        pairs.extend((i,int(j)) for j in js if np.isfinite(row[j]) and row[j]>=min_similarity)
    return np.asarray(pairs,dtype=int).reshape(-1,2)


def assignment(pa,pb,pairs,matrix,residual_px=6.):
    """Return unique proposal indices, maximizing count before residual quality."""
    if not len(pairs):return np.empty(0,int)
    sources,si=np.unique(pairs[:,0],return_inverse=True)
    targets,ti=np.unique(pairs[:,1],return_inverse=True)
    residual=np.linalg.norm(pa[pairs[:,0]]@matrix[:2,:2].T+matrix[:2,2]-pb[pairs[:,1]],axis=1)
    allowed=residual<residual_px
    # Any extra valid match outweighs the total possible valid residual cost.
    penalty=(len(sources)+1)*residual_px**2
    costs=np.full((len(sources),len(targets)+len(sources)),2*penalty)
    costs[:,len(targets):]=penalty
    indices=np.full((len(sources),len(targets)),-1,int)
    for index in np.flatnonzero(allowed):
        costs[si[index],ti[index]]=residual[index]**2
        indices[si[index],ti[index]]=index
    rr,cc=linear_sum_assignment(costs)
    legal=cc<len(targets)
    return indices[rr[legal],cc[legal]]


def two_point_hypotheses(pa,pb,pairs):
    """Vectorized closed-form orientation-preserving 2D similarity fits."""
    if len(pairs)<2:return np.empty((0,3,3)),0
    ij=np.asarray(list(combinations(range(len(pairs)),2)))
    first,second=pairs[ij[:,0]],pairs[ij[:,1]]
    da=pa[second[:,0]]-pa[first[:,0]];db=pb[second[:,1]]-pb[first[:,1]]
    legal=(first[:,0]!=second[:,0])&(first[:,1]!=second[:,1])
    denom=(da*da).sum(axis=1);legal&=(denom>1e-12)&((db*db).sum(axis=1)>1e-12)
    first=first[legal];da=da[legal];db=db[legal];denom=denom[legal]
    a=(da*db).sum(axis=1)/denom
    b=(da[:,0]*db[:,1]-da[:,1]*db[:,0])/denom
    matrix=np.repeat(np.eye(3)[None],len(a),axis=0)
    matrix[:,0,0]=a;matrix[:,0,1]=-b;matrix[:,1,0]=b;matrix[:,1,1]=a
    matrix[:,:2,2]=pb[first[:,1]]-np.einsum('nij,nj->ni',matrix[:,:2,:2],pa[first[:,0]])
    return matrix,len(ij)


def relational_filter(pa,pb,pairs,scores,residual_px=6.):
    pa=np.asarray(pa,float);pb=np.asarray(pb,float);pairs=np.asarray(pairs,int).reshape(-1,2)
    diagnostic=dict(status='insufficient_consensus',minimum_consensus=4,residual_threshold_histology_px=residual_px,
        proposals=len(pairs),hypotheses_tested=0,valid_pair_hypotheses=0,assignment_refinements=0)
    if len(pairs)<4:return None,np.empty((0,2),int),diagnostic
    matrices,attempts=two_point_hypotheses(pa,pb,pairs)
    diagnostic.update(hypotheses_tested=attempts,valid_pair_hypotheses=len(matrices))
    # Upper bound applies at INITIAL hypothesis, not after refitting. Pruning
    # is a heuristic for refined hypotheses; no global-optimality claim.
    # No source/target anatomical dimensions or frame centring used.
    upper=np.empty(len(matrices),int)
    for start in range(0,len(matrices),512):
        batch=matrices[start:start+512]
        pred=np.einsum('hij,pj->hpi',batch[:,:2,:2],pa[pairs[:,0]])+batch[:,None,:2,2]
        ok=np.linalg.norm(pred-pb[pairs[:,1]][None],axis=2)<residual_px
        for offset,legal in enumerate(ok):
            upper[start+offset]=min(len(np.unique(pairs[legal,0])),len(np.unique(pairs[legal,1])))
    best_key=(3,-np.inf,-np.inf);best=None;best_pairs=np.empty((0,2),int)
    for h in np.argsort(-upper,kind='stable'):
        if upper[h]<best_key[0]:break
        matrix=matrices[h]
        selected=assignment(pa,pb,pairs,matrix,residual_px)
        if len(selected)<4:continue
        chosen=pairs[selected]
        if not noncollinear(pa[chosen[:,0]],pb[chosen[:,1]]):continue
        diagnostic['assignment_refinements']+=1
        for _ in range(3):
            model=SimilarityTransform.from_estimate(pa[chosen[:,0]],pb[chosen[:,1]])
            if not isinstance(model,SimilarityTransform):break
            matrix=model.params
            new=pairs[assignment(pa,pb,pairs,matrix,residual_px)]
            if len(new)<4:chosen=new;break
            if np.array_equal(new,chosen):chosen=new;break
            chosen=new
        if len(chosen)<4 or not noncollinear(pa[chosen[:,0]],pb[chosen[:,1]]):continue
        residual=np.linalg.norm(pa[chosen[:,0]]@matrix[:2,:2].T+matrix[:2,2]-pb[chosen[:,1]],axis=1)
        # Last assignment is always checked against the returned matrix.
        key=(len(chosen),-float((residual**2).sum()),float(scores[chosen[:,0],chosen[:,1]].sum()))
        if key>best_key:best_key=key;best=matrix.copy();best_pairs=chosen.copy()
    if best is not None:diagnostic.update(status='estimated_model_not_anatomically_verified',retained=len(best_pairs))
    return best,best_pairs,diagnostic


def exhaustive_triples(pa,pb,pairs,residual_px=6.):
    """Diagnostic: enumerate original three-sample hypotheses, without truth.

    Same noncollinearity, min4 and final refit/residual checks as the baseline.
    Exhaustive sampling replaces RANSAC randomness; count/SSE selects hypotheses.
    """
    pairs=np.asarray(pairs,int);a=pa[pairs[:,0]];b=pb[pairs[:,1]]
    best=None;best_keep=np.zeros(len(pairs),bool);best_key=(3,-np.inf)
    valid=0;tested=0;largest_raw=0
    for indices in combinations(range(len(pairs)),3):
        tested+=1;indices=np.asarray(indices)
        if not noncollinear(a[indices],b[indices]):continue
        valid+=1;model=SimilarityTransform.from_estimate(a[indices],b[indices])
        if not isinstance(model,SimilarityTransform):continue
        residual=model.residuals(a,b);keep=residual<residual_px
        largest_raw=max(largest_raw,int(keep.sum()))
        if keep.sum()<4 or not noncollinear(a[keep],b[keep]):continue
        model=SimilarityTransform.from_estimate(a[keep],b[keep])
        if not isinstance(model,SimilarityTransform):continue
        residual=model.residuals(a,b);keep=residual<residual_px
        if keep.sum()<4 or not noncollinear(a[keep],b[keep]):continue
        key=(int(keep.sum()),-float((residual[keep]**2).sum()))
        if key>best_key:best_key=key;best=model.params.copy();best_keep=keep.copy()
    return best,pairs[best_keep],dict(triples_tested=tested,noncollinear_triples=valid,
        largest_raw_consensus=largest_raw,retained=int(best_keep.sum()),
        status='estimated_model_not_anatomically_verified' if best is not None else 'insufficient_consensus')
