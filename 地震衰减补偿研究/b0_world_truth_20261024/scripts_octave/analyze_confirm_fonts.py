"""Run the unchanged scripts/analyze_confirm.py with a CJK fallback font (the frozen font list names
Microsoft YaHei, which Linux lacks). Only the figure font changes."""
import runpy, sys
from pathlib import Path
import matplotlib
_update = matplotlib.RcParams.update
def _with_cjk(self, *a, **k):
    _update(self, *a, **k)
    fam = list(self['font.sans-serif'])
    if 'Noto Sans CJK SC' not in fam:
        dict.__setitem__(self, 'font.sans-serif', fam[:1] + ['Noto Sans CJK SC'] + fam[1:])
matplotlib.RcParams.update = _with_cjk
script = Path(__file__).resolve().parents[1] / 'scripts' / 'analyze_confirm.py'
sys.argv = [str(script)]
runpy.run_path(str(script), run_name='__main__')
