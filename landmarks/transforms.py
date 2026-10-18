"""2D similarity transforms in (column,row) coordinates, rows increasing down."""
import numpy as np
from scipy.ndimage import map_coordinates, gaussian_filter, minimum_filter


def apply_points(points_cr, matrix):
    p=np.asarray(points_cr,dtype=float)
    return (np.concatenate([p,np.ones((*p.shape[:-1],1))],axis=-1) @ matrix.T)[...,:2]


def similarity(shape_rc, clockwise_deg=20., scale=1.15, shift_cr=(40.,25.)):
    if not np.isfinite([clockwise_deg,scale,*shift_cr]).all() or scale<=0:
        raise ValueError('Transform parameters must be finite and scale positive')
    h,w=shape_rc
    centre=np.array([(w-1)/2,(h-1)/2])
    theta=np.deg2rad(clockwise_deg)
    # Positive angle is clockwise on a screen where rows increase downward.
    linear=scale*np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    matrix=np.eye(3); matrix[:2,:2]=linear
    matrix[:2,2]=centre-linear@centre+np.asarray(shift_cr)
    return matrix


def footprint(shape_rc):
    h,w=shape_rc
    return np.array([[-.5,-.5],[w-.5,-.5],[w-.5,h-.5],[-.5,h-.5]])


def expanded_transform(shape_rc, clockwise_deg=20., scale=1.15, shift_cr=(40.,25.), padding=10):
    if padding<0: raise ValueError('Padding must be nonnegative')
    base=similarity(shape_rc,clockwise_deg,scale,(0,0))
    requested=similarity(shape_rc,clockwise_deg,scale,shift_cr)
    # Preserve a reference frame: include original and unshifted rotated bounds,
    # with symmetric translation allowance. Do NOT tight-crop the shifted image.
    bounds=np.concatenate([footprint(shape_rc),apply_points(footprint(shape_rc),base)])
    margin=np.ceil(np.abs(shift_cr))+padding
    first=np.floor(bounds.min(axis=0)+.5)-margin
    last=np.ceil(bounds.max(axis=0)-.5)+margin
    offset=-first
    canvas=np.eye(3); canvas[:2,2]=offset
    effective=canvas@requested
    size=(last-first+1).astype(int)
    output_shape=(int(size[1]),int(size[0]))
    transformed=apply_points(footprint(shape_rc),effective)
    if np.any(transformed<-.5-1e-8) or np.any(transformed>size-.5+1e-8):
        raise ValueError('Transformed source footprint clipped')
    return effective,output_shape,offset,requested


def antialias_sigma(matrix):
    scales=np.linalg.svd(matrix[:2,:2],compute_uv=False)
    if not np.allclose(scales,scales[0],rtol=1e-7):
        raise ValueError('This sampler supports uniform similarity scale only')
    return float(.5*np.sqrt(max(0.,1/scales[0]**2-1)))


def warp(image,valid,matrix,output_shape,background=65535.):
    if image.shape!=valid.shape: raise ValueError('Mask shape differs from image')
    sigma=antialias_sigma(matrix)
    if sigma>1e-7:
        # Prefilter reductions; reject samples whose filter footprint touches
        # invalid support. Never let padding masquerade as source information.
        image=gaussian_filter(np.where(valid,image,background).astype(float),sigma=sigma,
                              mode='constant',cval=background,truncate=4.)
        radius=int(4*sigma+.5)
        valid=minimum_filter(valid.astype(np.uint8),size=2*radius+1,mode='constant',cval=0).astype(bool)
    rows,cols=np.indices(output_shape,dtype=float)
    source=apply_points(np.stack([cols,rows],axis=-1),np.linalg.inv(matrix))
    coords=np.stack([source[...,1],source[...,0]])
    output=map_coordinates(image,coords,order=1,mode='constant',cval=background,prefilter=False)
    # A sample is valid only if all its bilinear contributors have valid support.
    support=map_coordinates(valid.astype(float),coords,order=1,mode='constant',cval=0.,prefilter=False)
    inside=(source[...,0]>=0)&(source[...,0]<=image.shape[1]-1)&(source[...,1]>=0)&(source[...,1]<=image.shape[0]-1)
    output_valid=inside & (support>=1.-1e-12)
    output[~output_valid]=background
    return output,output_valid


def sample_parameters(seed,count):
    """Fixed initial experiment ranges; save sampled values as well as seed."""
    if not isinstance(count,int) or not 1<=count<=100:
        raise ValueError('Batch size must be between 1 and 100')
    if seed<0: raise ValueError('Seed must be nonnegative')
    rng=np.random.default_rng(seed)
    return [dict(angle=float(rng.uniform(-30,30)),scale=float(rng.uniform(.85,1.2)),
                 shift=list(map(int,rng.integers(-50,51,size=2)))) for _ in range(count)]
