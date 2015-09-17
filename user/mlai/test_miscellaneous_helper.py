from __future__ import division  #force floating point division
import unittest
import miscellaneous_helper as misc
import os

class TestMiscellaneousHelper(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_get_calling_function_name(self):
        function_name = misc.getCallingFunctionName()        
        self.assertEqual(function_name,'test_get_calling_function_name')

    def test_get_base_directory_of_calling_script(self):
        base_dir = misc.get_base_directory_of_calling_script()
        self.assertEqual(base_dir, 'utility')

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMiscellaneousHelper)
    unittest.TextTestRunner(verbosity=2).run(suite)