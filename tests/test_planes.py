import unittest
import numpy as np
from landmarks.planes import make_grid,extract,pixel_to_world,world_to_pixel

class PlaneTests(unittest.TestCase):
    def test_identity_slice_and_row_direction(self):
        data=np.arange(5*6*7,dtype=float).reshape(5,6,7)
        shape,m=make_grid([np.eye(4)],[data.shape],2.,1.)
        im,valid=extract(data,np.eye(4),shape,m)
        np.testing.assert_array_equal(im,data[:,2,:].T[::-1])
        self.assertTrue(valid.all())
        np.testing.assert_array_equal(pixel_to_world([[0,0],[4,6]],m),[[0,2,6],[4,2,0]])

    def test_oblique_source_sampling_and_pixel_roundtrip(self):
        a=np.array([[.7,.1,0,-3],[0,.8,.05,2],[.03,0,.9,-1],[0,0,0,1.]])
        ijk=np.moveaxis(np.indices((20,21,22),dtype=float),0,-1)
        world=ijk@a[:3,:3].T+a[:3,3]
        data=2*world[...,0]+3*world[...,1]-world[...,2]+7
        m=np.array([[.2,0,0],[0,0,8],[0,-.2,10],[0,0,1.]])
        im,valid=extract(data,a,(20,20),m)
        r,c=np.indices(im.shape)
        w=pixel_to_world(np.stack([c,r],-1),m)
        expected=2*w[...,0]+3*w[...,1]-w[...,2]+7
        self.assertTrue(valid.all())
        np.testing.assert_allclose(im,expected,atol=1e-11)
        np.testing.assert_allclose(world_to_pixel(w,m),np.stack([c,r],-1),atol=1e-12)
        with self.assertRaises(ValueError): world_to_pixel([[0,9,0]],m)

    def test_outside_support_is_invalid(self):
        m=np.array([[1,0,-1],[0,0,2],[0,-1,4],[0,0,1.]])
        im,valid=extract(np.ones((5,5,5)),np.eye(4),(5,7),m,background=99)
        self.assertFalse(valid[:,0].any()); self.assertFalse(valid[:,-1].any())
        np.testing.assert_array_equal(im[:,0],99)
        self.assertTrue(valid[:,1:-1].all())
        with self.assertRaises(ValueError): make_grid([np.eye(4)],[(5,5,5)],100)

if __name__=='__main__': unittest.main()
