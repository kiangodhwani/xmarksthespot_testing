"""Independent multiscale Harris candidates; no truth or counterpart input."""
import numpy as np
from scipy.ndimage import gaussian_filter, distance_transform_edt
from skimage.feature import corner_harris, corner_peaks


def detect_candidates(image, valid, tissue, scales=(2.,4.,8.), max_points=180,
                      min_distance=8., threshold_rel=.01, tiles=6, per_tile=10):
    image=np.asarray(image); valid=np.asarray(valid); tissue=np.asarray(tissue)
    if image.ndim!=2 or image.shape!=valid.shape or image.shape!=tissue.shape:
        raise ValueError('Image, valid and tissue arrays must share shape')
    values=image[valid]
    if not len(values) or not np.isfinite(values).all(): raise ValueError('Invalid intensities')
    lo,hi=np.percentile(values,[1,99])
    if hi<=lo: return np.empty((0,5)),dict(count=0,reason='No intensity variation')
    normal=np.clip((image.astype(float)-lo)/(hi-lo),0,1)
    # Fill unsupported values, but never multiply grayscale by tissue mask.
    normal[~valid]=float(np.median(normal[valid]))
    # Pad before distance transform so array edges count as unsupported too.
    safe_distance=distance_transform_edt(np.pad(valid,1))[1:-1,1:-1]
    allowed=valid & tissue
    candidates=[]; per_scale=[]
    for sigma in scales:
        smooth=gaussian_filter(normal,sigma=sigma/2,mode='nearest')
        response=corner_harris(smooth,method='k',k=.04,sigma=sigma)
        safe=allowed & (safe_distance>np.ceil(4*sigma+2*sigma+3))
        maximum=float(np.max(response[safe])) if safe.any() else 0.
        if maximum<=0: continue
        # Scores are relative per-image/per-scale responses, NOT probabilities.
        # Find genuine response maxima BEFORE filtering locations by the mask.
        # Zeroing outside tissue first could manufacture maxima at mask borders.
        response=np.maximum(response,0)/maximum
        coords=corner_peaks(response,min_distance=4,threshold_abs=threshold_rel,
                            threshold_rel=0,exclude_border=False)
        coords=coords[safe[coords[:,0],coords[:,1]]][:500]
        per_scale.append(dict(sigma_px=sigma,peaks=len(coords)))
        for r,c in coords: candidates.append((float(c),float(r),float(response[r,c]),float(sigma)))
    if not candidates: return np.empty((0,5)),dict(count=0,per_scale=per_scale)
    rr,cc=np.where(allowed); xmin,xmax=cc.min(),cc.max()+1; ymin,ymax=rr.min(),rr.max()+1
    occupancy={}; selected=[]
    for c,r,score,sigma in sorted(candidates,key=lambda p:(-p[2],p[1],p[0],p[3])):
        tile=(min(tiles-1,int((c-xmin)*tiles/(xmax-xmin))),min(tiles-1,int((r-ymin)*tiles/(ymax-ymin))))
        if occupancy.get(tile,0)>=per_tile: continue
        if any((c-p[0])**2+(r-p[1])**2<min_distance**2 for p in selected): continue
        selected.append((c,r,score,sigma,tile[1]*tiles+tile[0])); occupancy[tile]=occupancy.get(tile,0)+1
        if len(selected)>=max_points: break
    points=np.asarray(selected,dtype=float).reshape(-1,5)
    return points,dict(count=len(points),per_scale=per_scale,occupied_tiles=len(occupancy),
        scales_px=list(scales),max_points=max_points,min_spacing_px=min_distance,
        threshold_rel=threshold_rel,tile_grid=[tiles,tiles],per_tile_limit=per_tile,
        normalisation_percentiles=[float(lo),float(hi)],
        columns=['col_px','row_px','relative_response','sigma_px','tile_id'],
        border_policy='Exclude detector footprint near invalid support/image edges; tissue mask only filters candidate centres',
        limitation='Multiscale baseline; not intrinsically scale invariant, anatomically labelled, or cross-modality validated')
