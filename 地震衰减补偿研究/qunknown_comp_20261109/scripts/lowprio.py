"""Run a script at below-normal priority (Windows) with one BLAS thread: python lowprio.py script.py args...
Set OMP/OPENBLAS/MKL_NUM_THREADS=1 in the environment as well (they must be set before numpy is imported)."""
import sys, os, runpy, ctypes
if os.name == "nt":
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)   # BELOW_NORMAL
sys.argv = sys.argv[1:]; sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
runpy.run_path(sys.argv[0], run_name="__main__")
