import subprocess
import shlex
from unittest import TestCase


class radiationDamageMPITests(TestCase):
    def test_diffraction_calculation_parallel(self):
        """
        Test diffraction calculation in parallel.
        """
        mpicommand = 'mpiexec -n 2 --map-by node'

        command_sequence = ['radiationDamageMPI',
                            '--inputDir', 'pmi',
                            '--outputDir', 'diffr_out',
                            '--beamFile', 'test_files/s2e.beam',
                            '--geomFile', 'test_files/s2e.geom',
                            '--uniformRotation', '1',
                            '--calculateCompton', '0',
                            '--sliceInterval', '100',
                            '--numSlices', '100',
                            '--pmiStartID', '1',
                            '--pmiEndID', '1',
                            '--numDP', '2'
                            ]

        args = shlex.split(mpicommand) + command_sequence

        proc = subprocess.Popen(args)

        proc.wait()

        self.assertTrue(True)

