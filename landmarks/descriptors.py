"""Experimental unsigned polar-gradient descriptors; image and point inputs only.

Four radial rings x twelve angular sectors x eight unsigned orientation bins.
Orientation is measured relative to radial direction. Circular sector shifts
provide a discrete rotation search; multiple radii provide a discrete scale bank.
This is a classical baseline, not a claim of modality independence.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

RADII=(16.,24.,36.)


def describe(image,valid,points_cr,radii=RADII):
    image=np.asarray(image);valid=np.asarray(valid);points=np.asarray(points_cr,dtype=float)
    if image.ndim!=2 or image.shape!=valid.shape: raise ValueError('Image/support shape mismatch')
    values=image[valid]
    if not len(values) or not np.isfinite(values).all(): raise ValueError('Invalid support intensities')
    lo,hi=np.percentile(values,[1,99])
    normal=np.clip((image.astype(float)-lo)/max(hi-lo,1e-12),0,1)
    normal[~valid]=float(np.median(normal[valid]))
    side=41; half=20
    yy,xx=np.indices((side,side),dtype=float);xx-=half;yy-=half
    rr=np.hypot(xx,yy);angle=np.arctan2(yy,xx)
    disc=(rr>=1)&(rr<=18)
    ring=np.minimum((rr[disc]/18*4).astype(int),3)
    sector=np.mod((angle[disc]+np.pi)/(2*np.pi)*12,12).astype(int)
    weights=np.exp(-rr[disc]**2/(2*14**2))
    bank=np.zeros((len(points),len(radii),4,12,8),dtype=np.float32)
    usable=np.zeros((len(points),len(radii)),bool)
    for k,radius in enumerate(radii):
        if radius<=0: raise ValueError('Radius must be positive')
        # Prefilter before patch resampling; fixed rule independent of truth scale.
        sigma=max(.7,radius/half*.7)
        smooth=gaussian_filter(normal,sigma=sigma,mode='nearest')
        for i,(col,row) in enumerate(points):
            coords=np.stack([row+yy*radius/half,col+xx*radius/half])
            # Reject a patch if any contributing interpolation sample is invalid.
            support=map_coordinates(valid.astype(float),coords,order=1,mode='constant',cval=0,prefilter=False)
            if not (support>=1-1e-12).all(): continue
            patch=map_coordinates(smooth,coords,order=1,mode='nearest',prefilter=False)
            gy,gx=np.gradient(gaussian_filter(patch,.7))
            magnitude=np.hypot(gx,gy)[disc]*weights
            orientation=np.mod(np.arctan2(gy,gx)[disc]-angle[disc],np.pi)/np.pi*8
            bins=np.floor(orientation).astype(int)%8; fraction=orientation-np.floor(orientation)
            hist=np.zeros((4,12,8))
            np.add.at(hist,(ring,sector,bins),magnitude*(1-fraction))
            np.add.at(hist,(ring,sector,(bins+1)%8),magnitude*fraction)
            # Smooth angular discretisation; closed angular/orientation dimensions.
            hist=gaussian_filter(hist,sigma=(.35,.5,.5),mode=('nearest','wrap','wrap'))
            if hist.sum()<1e-8: continue
            hist=np.sqrt(hist/hist.sum())
            hist/=np.linalg.norm(hist)
            bank[i,k]=hist;usable[i,k]=True
    return bank,usable


def similarity(a,a_valid,b,b_valid):
    """Max cosine over independently computed radii and 12 cyclic sector shifts.

    Returns score matrix, NOT accepted correspondences. Invalid pairs get -1.
    """
    if a.shape[1:]!=b.shape[1:]: raise ValueError('Descriptor dimensions differ')
    n,scales=a.shape[:2];m=b.shape[0]
    scores=np.full((n,m),-1.,dtype=float)
    if n==0 or m==0: return scores
    left=a.reshape(n*scales,-1)
    legal=a_valid[:,:,None,None]&b_valid[None,None,:,:]
    for shift in range(a.shape[3]):
        right=np.roll(b,shift,axis=3).reshape(m*scales,-1)
        values=(left@right.T).reshape(n,scales,m,scales)
        values=np.where(legal,values,-1.)
        scores=np.maximum(scores,values.max(axis=(1,3)))
    return np.minimum(scores,1.)
