import unittest
import numpy as np
from landmarks.masking import otsu_mask

class MaskTests(unittest.TestCase):
    def test_padding_cannot_change_threshold_or_become_tissue(self):
        image=np.zeros((12,12)); image[3:9,3:9]=100
        valid=np.ones_like(image,bool)
        base,report=otsu_mask(image,valid,True,0)
        padded=np.pad(image,4,constant_values=1e12); support=np.pad(valid,4)
        result,new=otsu_mask(padded,support,True,0)
        self.assertEqual(report['threshold'],new['threshold'])
        np.testing.assert_array_equal(result['raw'][4:-4,4:-4],base['raw'])
        self.assertFalse((result['raw'] & ~support).any())

    def test_polarity_holes_specks_and_input_preservation(self):
        image=np.full((16,16),100.); image[3:12,3:12]=20
        image[6:8,6:8]=100  # cavity: must not fill
        image[14,14]=20  # isolated speck: only cleaned mask should remove it
        before=image.copy(); valid=np.ones_like(image,bool)
        result,_=otsu_mask(image,valid,False,1)
        self.assertTrue(result['raw'][14,14]); self.assertFalse(result['cleaned'][14,14])
        self.assertFalse(result['cleaned'][6:8,6:8].any())
        self.assertTrue(result['cleaned'][4,4])
        np.testing.assert_array_equal(image,before)
        with self.assertRaises(ValueError): otsu_mask(image,np.zeros_like(valid),False)


class HysteresisTests(unittest.TestCase):
    def test_pale_connection_cavity_and_detached_region(self):
        from landmarks.masking import histology_hysteresis
        image=np.full((30,30),100.)
        image[5:20,5:20]=99.  # pale connected tissue
        image[8:12,8:12]=20.  # strong seed
        image[14:16,14:16]=100.  # cavity, must remain excluded
        image[25:27,25:27]=99.  # isolated pale patch, not seeded
        original=image.copy(); valid=np.ones_like(image,bool)
        result,_=histology_hysteresis(image,valid)
        self.assertTrue(result['mask'][6,6])
        self.assertFalse(result['strong'][6,6])
        self.assertFalse(result['mask'][14:16,14:16].any())
        self.assertFalse(result['mask'][25:27,25:27].any())
        np.testing.assert_array_equal(image,original)

    def test_padding_does_not_change_result(self):
        from landmarks.masking import histology_hysteresis
        image=np.full((20,20),100.); image[5:15,5:15]=95.; image[8:12,8:12]=30.
        valid=np.ones_like(image,bool)
        base,report=histology_hysteresis(image,valid)
        padded=np.pad(image,5,constant_values=-1e9); support=np.pad(valid,5)
        result,new=histology_hysteresis(padded,support)
        self.assertEqual(report,new)
        np.testing.assert_array_equal(result['mask'][5:-5,5:-5],base['mask'])
        self.assertFalse((result['mask']&~support).any())

if __name__=='__main__': unittest.main()
