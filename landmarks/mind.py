"""2D MIND-style sparse context descriptor, not a reproduction of 3D MIND-SSC.

Self-similarity exp(-local patch SSD / local mean neighbour SSD), max-normalised.
Reference: Heinrich et al. 2012, doi:10.1016/j.media.2012.05.008.
All sampling settings are fixed independently of evaluation geometry.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, minimum_filter, map_coordinates

RADII = (16., 24., 36.)
ANGLES = tuple(range(0, 360, 30))
OFFSETS = ((4,0), (-4,0), (0,4), (0,-4), (4,4), (-4,4), (4,-4), (-4,-4))
CONFIG = dict(name='mind2d_context_v1', radii_working_px=RADII,
    angles_clockwise_degrees=ANGLES, canonical_side=41, context_centres=[-12,-6,0,6,12],
    offsets_col_row=OFFSETS, ssd_gaussian_sigma=1., ssd_kernel_side=5,
    presmoothing='sigma=max(0.7, radius/20*0.7), truncate=4; full support required',
    normalization='exp(-SSD/mean cardinal SSD), divide channel max, flatten and L2',
    similarity='max cosine: source angle 0, all target angles and all 3x3 radius pairs',
    caveat='2D MIND-style adaptation; not original dense 3D MIND or MIND-SSC; scores not probabilities')


def sample_patch(image, point, radius, angle):
    y,x=np.indices((41,41),dtype=float); x-=20; y-=20
    theta=np.deg2rad(angle); c,s=np.cos(theta),np.sin(theta)
    coords=np.array([point[1]+(s*x+c*y)*radius/20,
                     point[0]+(c*x-s*y)*radius/20])
    return map_coordinates(image,coords,order=1,mode='constant',cval=0,prefilter=False)


def patch_descriptor(patch):
    # Every retained context centre has full 5x5 and displaced-patch support.
    centres=np.array([8,14,20,26,32]); yy,xx=np.meshgrid(centres,centres,indexing='ij')
    distances=[]
    for dx,dy in OFFSETS:
        shifted=np.roll(patch,(-dy,-dx),axis=(0,1))
        ssd=gaussian_filter((patch-shifted)**2,1.,truncate=2.,mode='constant')
        distances.append(ssd[yy,xx])
    distances=np.stack(distances,axis=-1)
    variance=distances[...,:4].mean(axis=-1,keepdims=True)
    if variance.mean()<1e-10: return None
    variance=np.maximum(variance,max(1e-10,float(variance.mean())*.001))
    descriptor=np.exp(-distances/variance)
    descriptor/=np.maximum(descriptor.max(axis=-1,keepdims=True),1e-12)
    descriptor=descriptor.ravel()
    return descriptor/np.linalg.norm(descriptor)


def describe(image,valid,points_cr,radii=RADII):
    image=np.asarray(image,dtype=float); valid=np.asarray(valid,dtype=bool)
    points=np.asarray(points_cr,dtype=float).reshape(-1,2)
    if image.ndim!=2 or image.shape!=valid.shape or not valid.any():
        raise ValueError('Expected 2D image with nonempty matching support')
    if not np.isfinite(image[valid]).all() or not np.isfinite(points).all():
        raise ValueError('Nonfinite input')
    if not np.isfinite(radii).all() or np.any(np.asarray(radii)<=0):
        raise ValueError('Radii must be finite and positive')
    lo,hi=np.percentile(image[valid],[1,99])
    normal=np.clip((image-lo)/max(hi-lo,1e-12),0,1)
    normal[~valid]=0
    bank=np.zeros((len(points),len(radii),len(ANGLES),200),np.float32)
    usable=np.zeros(bank.shape[:3],bool)
    for k,radius in enumerate(radii):
        sigma=max(.7,radius/20*.7); reach=int(4*sigma+.5)
        smooth=gaussian_filter(normal,sigma,truncate=4.,mode='constant')
        safe=minimum_filter(valid.astype(np.uint8),size=2*reach+1,mode='constant',cval=0)
        for i,point in enumerate(points):
            for a,angle in enumerate(ANGLES):
                if not (sample_patch(safe.astype(float),point,radius,angle)>=1-1e-12).all(): continue
                descriptor=patch_descriptor(sample_patch(smooth,point,radius,angle))
                if descriptor is not None: bank[i,k,a]=descriptor;usable[i,k,a]=True
    return bank,usable


def similarity(a,av,b,bv,return_choices=False):
    n=len(a);m=len(b)
    scores=np.full((n,m),-1.,dtype=float)
    choices=np.full((n,m,3),-1,dtype=np.int16)
    for ka in range(a.shape[1]):
        for kb in range(b.shape[1]):
            for angle in range(b.shape[2]):
                values=a[:,ka,0]@b[:,kb,angle].T
                legal=av[:,ka,0,None]&bv[None,:,kb,angle]
                better=legal&(values>scores)
                scores[better]=values[better];choices[better]=[ka,kb,angle]
    scores=np.clip(scores,-1,1)
    return (scores,choices) if return_choices else scores
