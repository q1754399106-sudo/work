"""Raw MT19937 output words of MATLAB's rng(seed,'twister') (seed 0 -> 5489), as uint32 LE.
Usage: mlwords.py seed n outfile"""
import sys
import numpy as np
seed, n, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
seed = 5489 if seed == 0 else seed
np.random.RandomState(seed).randint(0, 2**32, size=n, dtype=np.uint32).astype("<u4").tofile(path)
