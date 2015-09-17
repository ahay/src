from __future__ import division  #force floating point division
import unittest
import numpy as np
import file_operation as file_op
import rsf.api


class TestFileOperation(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_read_rsf_header_for_given_dimension(self):
        rsf_file_name = "wz.02.H"
        rsf_input = rsf.api.Input(rsf_file_name)

        output = file_op.readRsfHeaderForGivenDimension(rsf_input,1)
        self.assertEqual(output["n1"],1075)
        self.assertAlmostEqual(output["d1"],0.004)
        self.assertEqual(output["o1"],0)
        self.assertEqual(output["label1"],"Time")
        self.assertEqual(output["unit1"],"s")

        output = file_op.readRsfHeaderForGivenDimension(rsf_input,2)
        self.assertEqual(output["n2"],120)
        self.assertAlmostEqual(output["d2"],1)
        self.assertEqual(output["o2"],0)
        self.assertEqual(output["label2"],"Trace")
        self.assertEqual(output["unit2"],None)
                
    def test_read_rsf_header_for_data_type_information(self):
        rsf_file_name = "wz.02.H"
        rsf_input = rsf.api.Input(rsf_file_name)

        output = file_op.readRsfHeaderForDataTypeInformation(rsf_input)
        self.assertEqual(output["in"],"stdin")
        self.assertEqual(output["esize"],4)
        self.assertEqual(output["type"],"float")
        self.assertEqual(output["form"],"xdr")


    def test_convert_rsf_header_into_dictionary(self):
        rsf_file_name = "wz.02.H"
        output = file_op.convertRsfHeaderIntoDictionary(rsf_file_name)
        
        self.assertEqual(output["n1"],1075)
        self.assertAlmostEqual(output["d1"],0.004)
        self.assertEqual(output["o1"],0)
        self.assertEqual(output["label1"],"Time")
        self.assertEqual(output["unit1"],"s")
        
        self.assertEqual(output["n2"],120)
        self.assertAlmostEqual(output["d2"],1)
        self.assertEqual(output["o2"],0)
        self.assertEqual(output["label2"],"Trace")
        self.assertEqual(output["unit2"],None)
        
        self.assertEqual(output["in"],"stdin")
        self.assertEqual(output["esize"],4)
        self.assertEqual(output["type"],"float")
        self.assertEqual(output["form"],"xdr")

    def test_get_array_dimensions_from_rsf_header(self):
        rsf_header = dict(n1=7, n2=3, n3=1)
        output = file_op.getArrayDimensionsFromRsfHeader(rsf_header)
        self.assertEqual(output, [7, 3, 1])

        rsf_header = dict(n1=3, n3=2)
        output = file_op.getArrayDimensionsFromRsfHeader(rsf_header)
        self.assertEqual(output, [3])

        rsf_header = dict(n1=1, n3=3, n2=2)
        output = file_op.getArrayDimensionsFromRsfHeader(rsf_header)
        self.assertEqual(output, [1, 2, 3])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFileOperation)
    unittest.TextTestRunner(verbosity=2).run(suite)