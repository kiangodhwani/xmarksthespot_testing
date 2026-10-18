import csv
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from export_all_landmarks import HEADER,landmark_rows,write_csv,compare_existing,score_metrics


class LandmarkExportTests(unittest.TestCase):
    def setUp(self):
        self.matrix=np.array([[.2,0,-10],[0,0,8],[0,-.2,30],[0,0,1]],float)

    def test_world_coordinates_row_col_and_prediction_preservation(self):
        pa=np.array([[0,0],[3,2]],float);pb=np.array([[100,200],[7,9]],float)
        pairs=np.array([[1,0],[0,1]])
        scores=np.array([[.3,.8],[.9,.2]])
        rows=landmark_rows(pa,pb,pairs,scores,self.matrix)
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[0][:2],[1,'candidate_1'])
        np.testing.assert_allclose(rows[0][2:],[-9.4,8,29.6,100,200,.9])
        # A distant/incorrect target is still exported: no evaluation-based filter.
        self.assertEqual(rows[1][5:],[7.,9.,.8])

    def test_header_only_failure_and_existing_csv_verification(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'landmarks.csv';write_csv(path,HEADER,[])
            with path.open(newline='') as stream:
                self.assertEqual(list(csv.reader(stream)),[HEADER])
            self.assertTrue(compare_existing(path,HEADER,[]))
            with self.assertRaises(ValueError):compare_existing(path,HEADER,[[1,'candidate_1',0,0,0,0,0,.8]])
        reference=dict(original_pixel_to_world_ras_mm=self.matrix.tolist(),
            transformed_pixel_to_world_ras_mm=self.matrix.tolist(),original_to_transformed=np.eye(3).tolist())
        metrics,errors=score_metrics(np.array([[0,0]]),np.array([[0,0]]),np.empty((0,2),int),reference,None,{})
        self.assertEqual(metrics['count'],0);self.assertEqual(len(errors),0)
        self.assertIsNone(metrics['fraction_within_2mm'])
        self.assertIsNone(metrics['random_candidate_pairing'])
        self.assertIsNone(metrics['transform_recovery'])

if __name__=='__main__':unittest.main()
