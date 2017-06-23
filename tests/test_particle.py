from unittest import TestCase
import numpy as np

from pysingfel.particle import *


class particleTests(TestCase):
    def test_readPDB(self):
        p = Particle()
        p.readPDB('test_files/1uf2.pdb', ff='WK')

        self.assertTrue(isinstance(p.ffTable, np.ndarray))

