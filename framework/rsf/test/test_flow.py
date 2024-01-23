import unittest
from rsf.flow import Flow, add_prefix

class TestFlow(unittest.TestCase):
    def test_add_prefix_with_prefix(self):
        output = add_prefix('spike', 'sf')
        self.assertEqual(output, 'sfspike')

    def test_add_prefix_without_prefix(self):
        output = add_prefix('sfspike', 'sf')
        self.assertEqual(output, 'sfspike')

    def test_without_source(self):
        output = Flow([], 'spike', '', stdin=0)
        self.assertEqual(output, 'sfspike > $TARGET')

    def test_with_source(self):
        output = Flow([], 'spike', '', stdin=1)
        self.assertEqual(output, '< $SOURCE sfspike > $TARGET')

    def test_null_target(self):
        output = Flow([], 'spike', '', stdin=0, stdout=0)
        self.assertEqual(output, 'sfspike >/dev/null')

    def test_program_with_prefix(self):
        output = Flow([], 'sfspike', '', stdin=0)
        self.assertEqual(output, 'sfspike > $TARGET')
    
    def test_mpi_program_without_source(self):
        output = Flow([], 'mpifwi', '', stdin=0)
        self.assertEqual(output, 'sfmpifwi --output=$TARGET')

    def test_mpi_program_with_source(self):
        output = Flow([], 'mpifwi', '', stdin=1)
        self.assertEqual(output, 'sfmpifwi --output=$TARGET --input=$SOURCE')
    
    def test_sfmpi(self):
        output = Flow([], 'mpi transp', '')
        self.assertEqual(output, 'sfmpi sftransp --output=$TARGET --input=$SOURCE')
    
    def test_bindir(self):
        output = Flow([], '/home/sfspike', '/home')
        self.assertEqual(output, '< $SOURCE /home/sfspike > $TARGET')

    def test_checkpar(self):
        with self.assertRaises(SystemExit) as cm:
            _ = Flow([], 'bandpass fli=10', '', checkpar=True)
        self.assertEqual(cm.exception.code, 1)
    
    def test_conjgrad_with_suffix(self):
        output = Flow([], 'conjgrad sfcausint', '', stdin=0)
        self.assertEqual(output, 'sfconjgrad sfcausint > $TARGET')

    def test_conjgrad_search_next_command(self):
        output = Flow([], 'conjgrad niter=10 sfcausint', '', stdin=0)
        self.assertEqual(output, 'sfconjgrad niter=10 sfcausint > $TARGET')

    def test_conjgrad_without_suffix(self):
        output = Flow([], 'conjgrad causint', '', stdin=0)
        self.assertEqual(output, 'sfconjgrad sfcausint > $TARGET')

    def test_conjgrad_with_mpi(self):
        output = Flow([], 'conjgradmpi mpifwi', '', stdin=0, mpirun='mpirun -np 16')
        self.assertEqual(output, 'sfconjgradmpi mpirun -np 16 sfmpifwi --output=$TARGET')
