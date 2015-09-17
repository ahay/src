from __future__ import division  #force floating point division
import unittest
import numpy as np
import math_operation as math_op

class TestMathOperation(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_find_lower_median_index_and_value(self):
                #even number of elements
        data = [3, 2, 7, 1]
        median_point = math_op.findLowerMedianIndexAndValue(data)
        self.assertEqual(median_point.index,1)
        self.assertEqual(median_point.value,2)
                
                #odd number of elements
        data = [7, 20, 8]
        median_point = math_op.findLowerMedianIndexAndValue(data)
        self.assertEqual(median_point.index,2)
        self.assertEqual(median_point.value,8)
     
                #negative elements
        data = [-1, -19, -3, 2.3]
        median_point = math_op.findLowerMedianIndexAndValue(data)
        self.assertEqual(median_point.index,2)
        self.assertEqual(median_point.value,-3)	
                

    def test_find_lower_median_trace_index_time_index_value_of_family_of_trace(self):    
                #even number of traces
        data_by_trace_time = np.array([[1, 3, 8, 2], [9, 12, 11, 2.1]])
        family_trace_median = math_op.findLowerMedianTraceIndexTimeIndexValueOfFamilyOfTrace(data_by_trace_time)
        self.assertEqual(family_trace_median.time_index,1)
        self.assertEqual(family_trace_median.trace_index,0)
        self.assertEqual(family_trace_median.value,3)
 
                #odd number of traces
        data_by_trace_time = np.array([[8.5, 8, 5], 
                                   [1, 9, 2],
                                   [10, 5.5, 7]])
        family_trace_median = math_op.findLowerMedianTraceIndexTimeIndexValueOfFamilyOfTrace(data_by_trace_time)
        self.assertEqual(family_trace_median.time_index,2)
        self.assertEqual(family_trace_median.trace_index,2)
        self.assertEqual(family_trace_median.value,7)
        
    def test_calculate_semblance(self):
        data = [1, 2, -3, 4]
        semblance = math_op.calculateSemblance(data)
        self.assertEqual(semblance,16/(4*30))
         
                #constant vector
        random_int = np.random.randint(1,100)
        data = random_int * np.ones([10, 1])
        semblance = math_op.calculateSemblance(data)
        self.assertEqual(semblance,1)
    
    def test_weight_signal_by_domain_exponentiated(self):
        list_domain = np.array([1, 2, 3])
        domain_exponent_power = 2
        value_by_trace_time = np.array([[4,5,6],[7,8,9]])
        out = math_op.weightSignalByDomainExponentiated(list_domain, value_by_trace_time, domain_exponent_power)
        expected_out = np.array([[4,20,54],[7,32,81]])
        np.testing.assert_array_equal(expected_out, out)

        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMathOperation)
    unittest.TextTestRunner(verbosity=2).run(suite)