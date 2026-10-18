import unittest
import numpy as np
from landmarks.descriptors import describe,similarity

class DescriptorTests(unittest.TestCase):
    def test_identity_contrast_reversal_and_geometry(self):
        y,x=np.indices((121,121),dtype=float)
        image=np.exp(-((x-49)**2+(y-57)**2)/55)+.8*np.exp(-((x-72)**2+(y-69)**2)/25)
        image+=.3*((x>60)&(y<52))
        points=np.array([[60.,60.]])
        valid=np.ones_like(image,bool)
        original=image.copy()
        a,av=describe(image,valid,points,radii=(24.,))
        b,bv=describe(1-image,valid,points,radii=(24.,))
        c,cv=describe(np.rot90(image),valid,points,radii=(24.,))
        self.assertTrue(av.all())
        self.assertAlmostEqual(similarity(a,av,a,av)[0,0],1.,places=5)
        self.assertGreater(similarity(a,av,b,bv)[0,0],.999)
        self.assertGreater(similarity(a,av,c,cv)[0,0],.95)
        np.testing.assert_array_equal(image,original)

    def test_invalid_patch_rejected_and_constant_has_no_signature(self):
        image=np.zeros((101,101));image[40:60,40:60]=1
        valid=np.ones_like(image,bool)
        d,v=describe(image,valid,[[2.,2.]],radii=(16.,))
        self.assertFalse(v.any())
        self.assertEqual(similarity(d,v,d,v)[0,0],-1)
        d,v=describe(np.ones_like(image),valid,[[50.,50.]],radii=(16.,))
        self.assertFalse(v.any())

if __name__=='__main__': unittest.main()
