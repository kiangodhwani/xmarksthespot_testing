import unittest
import numpy as np
from landmarks.candidates import detect_candidates
from landmarks.repeatability import repeatability

class CandidateTests(unittest.TestCase):
    def test_candidate_coordinates_and_no_mask_edge_detection(self):
        image=np.zeros((160,160)); image[40:90,45:105]=1
        valid=np.ones_like(image,bool); tissue=np.ones_like(valid)
        points,_=detect_candidates(image,valid,tissue,scales=(2.,),min_distance=8)
        self.assertGreaterEqual(len(points),4)
        for p in points:
            self.assertTrue(tissue[int(p[1]),int(p[0])])
        # A tissue-mask edge alone on a constant image must not create corners.
        mask=np.zeros_like(valid); mask[40:90,45:105]=True
        got,_=detect_candidates(np.ones_like(image),valid,mask)
        self.assertEqual(len(got),0)

    def test_invalid_padding_independence(self):
        image=np.zeros((160,160)); image[50:100,50:100]=1
        valid=np.ones_like(image,bool); tissue=valid.copy()
        base,_=detect_candidates(image,valid,tissue,scales=(2.,))
        got,_=detect_candidates(np.pad(image,20,constant_values=1e9),np.pad(valid,20),np.pad(tissue,20),scales=(2.,))
        np.testing.assert_allclose(got[:,:2]-20,base[:,:2])

    def test_mask_cannot_manufacture_response_peaks(self):
        image=np.zeros((160,160)); image[40:100,40:100]=1
        valid=np.ones_like(image,bool)
        # No true square corner lies inside this inset search mask.
        tissue=np.zeros_like(valid); tissue[45:95,45:95]=True
        points,_=detect_candidates(image,valid,tissue,scales=(2.,))
        self.assertEqual(len(points),0)

    def test_one_to_one_and_tolerance(self):
        report=repeatability(np.array([[0.,0.],[.1,0.]]),np.array([[0.,0.]]),(.2,))[0]
        self.assertEqual(report['one_to_one_pairs'],1)
        self.assertEqual(report['fraction_of_source'],.5)
        report=repeatability(np.array([[0.,0.]]),np.array([[3.,4.]]),(2.,5.))
        self.assertEqual([r['one_to_one_pairs'] for r in report],[0,1])

if __name__=='__main__': unittest.main()
