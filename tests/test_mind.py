import unittest
import numpy as np
from scipy.ndimage import map_coordinates
from landmarks.mind import describe, similarity


class MindTests(unittest.TestCase):
    def setUp(self):
        y,x=np.indices((151,151),dtype=float)
        self.image=np.exp(-((x-64)**2+(y-72)**2)/55)+.8*np.exp(-((x-85)**2+(y-84)**2)/25)
        self.image+=.3*((x>75)&(y<67))
        self.valid=np.ones_like(self.image,bool)
        self.point=[[75.,75.]]

    def test_identity_contrast_and_rotation_search(self):
        original=self.image.copy()
        a,av=describe(self.image,self.valid,self.point,radii=(24.,))
        b,bv=describe(1-self.image,self.valid,self.point,radii=(24.,))
        c,cv=describe(np.rot90(self.image),self.valid,self.point,radii=(24.,))
        self.assertTrue(av.all())
        self.assertAlmostEqual(similarity(a,av,a,av)[0,0],1,places=5)
        self.assertGreater(similarity(a,av,b,bv)[0,0],.999)
        self.assertGreater(similarity(a,av,c,cv)[0,0],.999)
        np.testing.assert_array_equal(self.image,original)

    def test_scale_search(self):
        y,x=np.indices(self.image.shape,dtype=float)
        enlarged=map_coordinates(self.image,[(y-75)/1.5+75,(x-75)/1.5+75],order=1)
        a,av=describe(self.image,self.valid,self.point,radii=(24.,))
        b,bv=describe(enlarged,self.valid,self.point,radii=(36.,))
        self.assertGreater(similarity(a,av,b,bv)[0,0],.98)

    def test_invalid_constant_and_padding_independence(self):
        d,v=describe(np.ones_like(self.image),self.valid,self.point,radii=(24.,))
        self.assertFalse(v.any())
        self.assertEqual(similarity(d,v,d,v)[0,0],-1)
        d,v=describe(self.image,self.valid,[[1.,1.]],radii=(24.,))
        self.assertFalse(v.any())
        valid=self.valid.copy();valid[:20]=False
        modified=self.image.copy();modified[~valid]=1e20
        a,av=describe(self.image,valid,self.point,radii=(24.,))
        b,bv=describe(modified,valid,self.point,radii=(24.,))
        np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(av,bv)
        # Gaussian prefilter support, not just patch centres, must be valid.
        valid[49,75]=False
        _,flags=describe(self.image,valid,self.point,radii=(24.,))
        self.assertFalse(flags[0,0,0])

    def test_empty_bank_and_choice_reconstructs_score(self):
        a,av=describe(self.image,self.valid,self.point,radii=(16.,24.))
        scores,choices=similarity(a,av,a,av,return_choices=True)
        ka,kb,angle=choices[0,0]
        self.assertAlmostEqual(scores[0,0],float(a[0,ka,0]@a[0,kb,angle]),places=5)
        empty,ev=describe(self.image,self.valid,[],radii=(16.,24.))
        self.assertEqual(similarity(empty,ev,a,av).shape,(0,1))


if __name__=='__main__': unittest.main()
