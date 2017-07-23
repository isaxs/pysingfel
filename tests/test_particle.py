import unittest
import numpy as np

from pysingfel.particle import *


class particleTests(unittest.TestCase):
    def test_readPDB(self):
        p = Particle()
        p.readPDB('test_files/5g3x.pdb', ff='WK')

        self.assertTrue(isinstance(p.ffTable, np.ndarray))

if __name__ == '__main__':
    unittest.main()
