import unittest
from rsf.flow import Flow

class TestFlow(unittest.TestCase):
    def test_flow(self):
        output = Flow([], 'spike', '', stdin=0)
        self.assertEqual(output, 'sfspike > $TARGET')