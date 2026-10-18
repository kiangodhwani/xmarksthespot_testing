"""Coronal sampling. Images index [row,col]; points use (col,row)."""
import numpy as np
from scipy.ndimage import map_coordinates


def make_grid(affines, shapes, y_mm, spacing_mm=0.2):
    if not np.isfinite([y_mm, spacing_mm]).all() or spacing_mm <= 0:
        raise ValueError('Plane and positive spacing must be finite')
    corners = []
    for affine, shape in zip(affines, shapes):
        if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3,:3])) < 1e-12:
            raise ValueError('Invalid source affine')
        if not np.allclose(affine[:3,:3], np.diag(np.diag(affine[:3,:3])), atol=1e-7):
            raise ValueError('This initial grid builder requires axis-aligned source volumes')
        indices = np.array(np.meshgrid(*[(0,n-1) for n in shape], indexing='ij')).reshape(3,-1).T
        world = indices @ affine[:3,:3].T + affine[:3,3]
        if not world[:,1].min() <= y_mm <= world[:,1].max():
            raise ValueError('Chosen plane is outside at least one source volume')
        corners.append(world)
    corners=np.concatenate(corners)
    lo=np.floor(corners.min(axis=0)/spacing_mm)*spacing_mm
    hi=np.ceil(corners.max(axis=0)/spacing_mm)*spacing_mm
    width=int(round((hi[0]-lo[0])/spacing_mm))+1
    height=int(round((hi[2]-lo[2])/spacing_mm))+1
    # [col,row,1] -> [x,y,z,1], exact plane y=y_mm.
    matrix=np.array([[spacing_mm,0,lo[0]], [0,0,y_mm],
                     [0,-spacing_mm,hi[2]], [0,0,1.]],dtype=float)
    return (height,width), matrix


def pixel_to_world(points_cr, matrix):
    points=np.asarray(points_cr,dtype=float)
    homogeneous=np.concatenate([points,np.ones((*points.shape[:-1],1))],axis=-1)
    return (homogeneous @ matrix.T)[...,:3]


def world_to_pixel(world_xyz, matrix, tolerance=1e-7):
    world=np.asarray(world_xyz,dtype=float)
    if not np.allclose(world[...,1],matrix[1,2],atol=tolerance,rtol=0):
        raise ValueError('Point is not on the saved coronal plane')
    return np.stack([(world[...,0]-matrix[0,2])/matrix[0,0],
                     (world[...,2]-matrix[2,2])/matrix[2,1]],axis=-1)


def extract(data, affine, shape_rc, matrix, background=0.):
    rows,cols=np.indices(shape_rc,dtype=float)
    world=pixel_to_world(np.stack([cols,rows],axis=-1),matrix)
    inverse=np.linalg.inv(affine)
    indices=world @ inverse[:3,:3].T + inverse[:3,3]
    valid=np.all((indices >= 0) & (indices <= np.asarray(data.shape)-1),axis=-1)
    image=map_coordinates(data,indices.reshape(-1,3).T,order=1,
                          mode='constant',cval=float(background),prefilter=False).reshape(shape_rc)
    if not valid.any():
        raise ValueError('Plane has no valid samples')
    return image,valid
