import unittest
import numpy as np
from landmarks.matching import propose,geometric_filter
from landmarks.transforms import similarity,apply_points

class MatchingTests(unittest.TestCase):
    def test_mutual_distinctiveness_and_empty_handling(self):
        mutual,accepted,ratios=propose(np.array([[.98,.3,.2],[.2,.97,.1],[.2,.3,.96]]))
        np.testing.assert_array_equal(accepted,[[0,0],[1,1],[2,2]])
        self.assertTrue((ratios<.9).all())
        _,accepted,_=propose(np.full((3,3),.9))
        self.assertEqual(len(accepted),0)
        self.assertEqual(len(propose(np.zeros((0,3)))[1]),0)
        self.assertEqual(len(propose(np.full((3,3),-1.))[1]),0)

    def test_geometric_recovery_rejects_outliers(self):
        rng=np.random.default_rng(9);source=rng.uniform(20,300,(25,2))
        matrix=similarity((350,350),17,1.1,(20,-7))
        target=apply_points(source,matrix)
        target[-7:]=rng.uniform(0,350,(7,2))
        fitted,keep,report=geometric_filter(source,target,residual_px=.1)
        self.assertIsNotNone(fitted);self.assertTrue(keep[:18].all());self.assertFalse(keep[18:].any())
        np.testing.assert_allclose(fitted,matrix,atol=1e-8)
        a,k,_=geometric_filter(source[:3],target[:3]);self.assertIsNone(a);self.assertFalse(k.any())

    def test_degenerate_input_not_accepted(self):
        source=np.stack([np.arange(10),np.zeros(10)],axis=1)
        model,keep,report=geometric_filter(source,source+4)
        self.assertIsNone(model);self.assertFalse(keep.any())

if __name__=='__main__': unittest.main()
