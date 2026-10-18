import unittest
import numpy as np
from landmarks.transforms import similarity,expanded_transform,apply_points,footprint,warp

class TransformTests(unittest.TestCase):
    def test_identity_and_asymmetric_translation(self):
        image=np.zeros((7,9)); image[1,2]=3; image[4,6]=7
        valid=np.ones_like(image,dtype=bool)
        identity=similarity(image.shape,0,1,(0,0))
        got,mask=warp(image,valid,identity,image.shape,background=0)
        np.testing.assert_array_equal(got,image); self.assertTrue(mask.all())
        matrix=similarity(image.shape,0,1,(3,2))
        got,mask=warp(image,valid,matrix,(12,15),background=0)
        self.assertEqual(got[3,5],3); self.assertEqual(got[6,9],7)
        self.assertFalse(mask[0,0]); self.assertEqual(np.count_nonzero(got),2)

    def test_clockwise_rotation_direction(self):
        image=np.zeros((7,7)); image[3,5]=9
        matrix=similarity(image.shape,90,1,(0,0))
        # Right of centre moves below centre for clockwise rotation.
        np.testing.assert_allclose(apply_points([[5,3]],matrix),[[3,5]],atol=1e-12)
        got,_=warp(image,np.ones_like(image,dtype=bool),matrix,image.shape,0)
        self.assertAlmostEqual(got[5,3],9)

    def test_canvas_and_roundtrips(self):
        points=np.array([[0,0],[699,630],[320.4,210.8]])
        for angle,scale,shift in [(0,1,(0,0)),(0,1,(40,25)),(20,1,(0,0)),(0,1.15,(0,0)),(20,1.15,(40,25)),(-37,.8,(-47,31))]:
            t,shape,offset,requested=expanded_transform((631,700),angle,scale,shift)
            p=apply_points(points,t)
            np.testing.assert_allclose(apply_points(p,np.linalg.inv(t)),points,atol=1e-10)
            np.testing.assert_allclose(p,apply_points(points,requested)+offset,atol=1e-10)
            corners=apply_points(footprint((631,700)),t)
            self.assertTrue((corners>=-.5).all())
            self.assertTrue((corners<=[shape[1]-.5,shape[0]-.5]).all())
        with self.assertRaises(ValueError): similarity((7,7),scale=0)

    def test_invalid_interpolation_contributor_is_excluded(self):
        image=np.ones((5,5)); valid=np.ones((5,5),bool); valid[2,2]=False
        t=similarity(image.shape,0,1,(.5,0))
        got,mask=warp(image,valid,t,image.shape,99)
        self.assertFalse(mask[2,2]); self.assertFalse(mask[2,3])
        self.assertEqual(got[2,3],99)



class RandomPairTests(unittest.TestCase):
    def test_seed_bounds_and_random_geometry(self):
        from landmarks.transforms import sample_parameters
        a=sample_parameters(42,10)
        self.assertEqual(a,sample_parameters(42,10))
        self.assertNotEqual(a,sample_parameters(43,10))
        for item in a:
            self.assertTrue(-30<=item['angle']<=30)
            self.assertTrue(.85<=item['scale']<=1.2)
            self.assertTrue(all(-50<=x<=50 for x in item['shift']))
            t,shape,offset,requested=expanded_transform((41,57),clockwise_deg=item['angle'],scale=item['scale'],shift_cr=item['shift'])
            p=footprint((41,57));q=apply_points(p,t)
            np.testing.assert_allclose(apply_points(q,np.linalg.inv(t)),p,atol=1e-10)
            self.assertTrue((q>=-.5).all());self.assertTrue((q<=[shape[1]-.5,shape[0]-.5]).all())

    def test_downscale_prefilter_and_input_preservation(self):
        from landmarks.transforms import antialias_sigma
        image=(np.indices((41,41)).sum(axis=0)%2).astype(float)
        original=image.copy();valid=np.ones_like(image,bool)
        t=similarity(image.shape,0,.85,(0,0))
        result,support=warp(image,valid,t,image.shape,0)
        self.assertGreater(antialias_sigma(t),0)
        self.assertEqual(antialias_sigma(np.eye(3)),0)
        self.assertTrue(np.isfinite(result).all());self.assertTrue(support.any())
        np.testing.assert_array_equal(image,original)

if __name__=='__main__': unittest.main()
