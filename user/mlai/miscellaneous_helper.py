'''
Created on Mar 26, 2015

@author: mlai
'''
from __future__ import division  #force floating point division
import numpy as np
import math
import collections
import logging
import sys
import os
import time
import inspect
import unittest
import test_miscellaneous_helper


def getCallingFunctionName():
    return inspect.stack()[1][3]

def createUniqueToFunctionCallLogger():
    current_time = time.strftime("%d%m%Y_%H%M%S") 
    rand_int_str = str(np.random.randint(0,10000))               
    unique_to_function_call_logger = logging.getLogger(current_time + rand_int_str) 
    unique_to_function_call_logger.setLevel(logging.DEBUG) 
    return unique_to_function_call_logger

def setupHandlerAndAddToLogger(handler, logger, logging_level = logging.DEBUG):
    handler.setLevel(logging_level)
    
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)

def mergeDictionaries(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result    

def get_directory_of_calling_script():
    return os.path.dirname(os.path.abspath(inspect.stack()[2][1]))

def get_base_directory_of_calling_script():
    dir_of_calling_script = get_directory_of_calling_script()
    head, tail = os.path.split(dir_of_calling_script)
    return tail

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(test_miscellaneous_helper.TestMiscellaneousHelper)
    unittest.TextTestRunner(verbosity=2).run(suite)