import unittest
import numpy as np
from landmarks.matching import propose
from landmarks.proposal_ablation import select_variants


class ProposalAblationTests(unittest.TestCase):
    def test_removes_each_condition_independently(self):
        # First mutual pair is high-score but ambiguous; second is distinctive
        # but below score threshold. The third passes both conditions.
        scores=np.full((6,6),-1.)
        scores[:2,:2]=[[.9,.89],[.2,.1]]
        scores[2:4,2:4]=[[.7,.1],[.05,.0]]
        scores[4:6,4:6]=[[.95,.1],[.05,.0]]
        original=scores.copy()
        mutual,selected,audit=select_variants(scores)
        np.testing.assert_array_equal(mutual,[[0,0],[2,2],[4,4]])
        np.testing.assert_array_equal(selected['current'],[[4,4]])
        np.testing.assert_array_equal(selected['similarity_only'],[[0,0],[4,4]])
        np.testing.assert_array_equal(selected['all_mutual'],mutual)
        np.testing.assert_array_equal(selected['current'],propose(scores)[1])
        self.assertFalse(audit[0]['passes_ratio'])
        self.assertFalse(audit[1]['passes_similarity'])
        np.testing.assert_array_equal(scores,original)

    def test_ties_invalid_and_single_alternative(self):
        for scores in [np.empty((0,3)),np.full((2,2),-1.)]:
            mutual,selected,audit=select_variants(scores)
            self.assertEqual(mutual.shape,(0,2));self.assertEqual(audit,[])
            for pairs in selected.values():self.assertEqual(pairs.shape,(0,2))
        _,selected,audit=select_variants(np.ones((1,1)))
        self.assertEqual(len(selected['current']),0)
        self.assertEqual(len(selected['all_mutual']),1)
        self.assertFalse(audit[0]['ratio_defined'])
        _,selected,audit=select_variants(np.ones((2,2)))
        self.assertEqual(len(selected['current']),0)
        self.assertFalse(audit[0]['ratio_defined'])

if __name__=='__main__':unittest.main()
