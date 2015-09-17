from __future__ import division  #force floating point division
import numpy as np
import collections
import scipy as Sci
import scipy.linalg
import unittest
import test_math_operation  

# This class provides the functionality we want. You only need to look at
# this if you want to know how this works. It only needs to be defined
# once, no need to muck around with its internals.
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

def divideEachRowByItsMaxAbsValue(input_by_row_column):
    assert input_by_row_column.ndim == 2
    max_of_each_row = np.max(np.abs(input_by_row_column),axis=1)
    if type(max_of_each_row) not in (np.ndarray,):
        max_of_each_row = np.array([max_of_each_row])
                #divide by eps if max of row is zero
    max_of_each_row[max_of_each_row == 0] = np.spacing(1)
    return input_by_row_column / max_of_each_row
 
def calculateLinearlySpacedValues(start_value, step_size, num_value):
    stop_value = start_value + step_size*num_value 
    return np.arange(start_value, stop_value ,step_size) 

# Based on a post by Peter Otten on comp.lang.python 2004/09/04.
# This uses the 'key' option to sort and the 'sorted' function, both
# new in Python 2.4.
def permutation_indices(data):
    return sorted(range(len(data)), key = data.__getitem__)


def findLowerMedianIndexAndValue(data):
    '''Suppose input has N elements, indexed by 0,1,2,...,N-1. 
       If N is odd, then  function returns (N-1)/2 element of input as value 
       and (N-1)/2 as index.
       If N is even, then function returns (N/2-1) element of input as value 
       and (N/2-1) as index. This is the "lower" median value/index pair.
    '''
    num_value = len(data)
    indicies_sorted = permutation_indices(data)
    point = collections.namedtuple('MedianPoint', ['index', 'value'])
    if (num_value % 2 == 0): #even 
        half_num_value = num_value/2 - 1        
    else: #odd
        half_num_value = (num_value - 1)/ 2
    
    point.index = indicies_sorted[int(half_num_value)]  
    point.value = data[point.index] 
    return point
        
def findLowerMedianTraceIndexTimeIndexValueOfFamilyOfTrace(data_by_trace_time): 
    '''finds the median value of family of traces along with the trace index and 
       and time index where median occurred
    '''
    assert data_by_trace_time.ndim == 2

    num_trace = data_by_trace_time.shape[0]
    num_time = data_by_trace_time.shape[1]
    
            #ravel returns a view of the original array whenever possible
    data_as_vector = np.ravel(data_by_trace_time)
    output = collections.namedtuple('FamilyTraceMedianPoint', ['trace_index','time_index', 'value'])

            #find median and location
    point = findLowerMedianIndexAndValue(data_as_vector)
    output.value = point.value
    median_coordinate_in_data_by_trace_time = np.unravel_index(point.index,
                                                                (num_trace,
                                                                 num_time))
    output.trace_index = median_coordinate_in_data_by_trace_time[0]
    output.time_index = median_coordinate_in_data_by_trace_time[1]
    return output

def calculateValuePercentile(start_value, end_value, value):
    spread_of_values = end_value - start_value
    return (value - start_value) / (spread_of_values)    

def calculateSemblance(data_vector):
    '''measures inputs correlation with a constant vector'''
    num_entry = np.size(data_vector)
    numerator = np.square(np.sum(data_vector))
    denominator = num_entry * np.sum(np.square(data_vector))
    return numerator / denominator

def weightSignalByDomainExponentiated(list_domain,
                                      value_by_trace_time,
                                      domain_exponent_power):

    list_domain_exponentiated = np.power(list_domain,domain_exponent_power)
    return np.multiply( 
           list_domain_exponentiated,
           value_by_trace_time)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(test_math_operation.TestMathOperation)
    unittest.TextTestRunner(verbosity=2).run(suite)
    