import unittest
import numpy as np
from landmarks.relational import topk_proposals,assignment,relational_filter,exhaustive_triples


class RelationalTests(unittest.TestCase):
    def test_topk_stable_invalid_and_threshold(self):
        s=np.array([[.8,.8,.75,.7],[-1,-1,-1,-1]])
        np.testing.assert_array_equal(topk_proposals(s),[[0,0],[0,1],[0,2]])
        self.assertEqual(topk_proposals(np.empty((0,4))).shape,(0,2))

    def test_alternatives_recover_transform_without_duplicate_votes(self):
        pa=np.array([[0,0],[30,0],[0,40],[50,60],[80,20]],float)
        theta=.3;r=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])*1.1
        pb=pa@r.T+[20,17];pb=np.vstack([pb,[[200,200],[230,190],[250,280],[190,300],[260,250]]])
        scores=np.zeros((5,10));scores[np.arange(5),np.arange(5)]=.9
        scores[np.arange(5),np.arange(5)+5]=.95
        proposals=topk_proposals(scores)
        model,chosen,diag=relational_filter(pa,pb,proposals,scores)
        self.assertIsNotNone(model)
        self.assertEqual(len(chosen),5)
        self.assertEqual(len(np.unique(chosen[:,0])),5);self.assertEqual(len(np.unique(chosen[:,1])),5)
        self.assertLess(np.max(np.linalg.norm(pa@model[:2,:2].T+model[:2,2]-pb[:5],axis=1)),1e-8)

    def test_assignment_maximizes_cardinality_and_allows_unmatched(self):
        pa=np.array([[0,0],[1,0],[100,100]],float);pb=np.array([[0,0],[5,0]],float)
        pairs=np.array([[0,0],[0,1],[1,0],[2,0]])
        chosen=pairs[assignment(pa,pb,pairs,np.eye(3))]
        np.testing.assert_array_equal(chosen,[[0,1],[1,0]])

    def test_exhaustive_known_model_and_degenerate_rejection(self):
        pa=np.array([[0,0],[20,0],[0,30],[50,60],[90,30]],float);pb=pa+[15,20]
        pairs=np.column_stack([np.arange(5),np.arange(5)])
        model,chosen,diag=exhaustive_triples(pa,pb,pairs)
        self.assertEqual(diag['triples_tested'],10);self.assertEqual(len(chosen),5)
        self.assertIsNotNone(model)
        line=np.column_stack([np.arange(5)*20,np.zeros(5)])
        self.assertIsNone(relational_filter(line,line,pairs,np.eye(5))[0])
        self.assertIsNone(exhaustive_triples(line,line,pairs)[0])

if __name__=='__main__':unittest.main()
