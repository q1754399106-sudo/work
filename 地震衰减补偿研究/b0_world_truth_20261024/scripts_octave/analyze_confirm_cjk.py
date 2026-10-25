"""Run the unchanged scripts/analyze_confirm.py with a CJK fallback font for the figure (the frozen
font list names Microsoft YaHei, which Linux lacks). Replaces analyze_confirm_fonts.py, whose patch
raised KeyError on matplotlib's partial RcParams objects before any scoring ran."""
import runpy, sys
from pathlib import Path
import matplotlib
_update = matplotlib.RcParams.update


def _with_cjk(self, *a, **k):
    _update(self, *a, **k)
    if 'font.sans-serif' in self.keys():
        fam = list(dict.__getitem__(self, 'font.sans-serif'))
        if 'Noto Sans CJK SC' not in fam:
            dict.__setitem__(self, 'font.sans-serif', fam[:1] + ['Noto Sans CJK SC'] + fam[1:])


matplotlib.RcParams.update = _with_cjk
script = Path(__file__).resolve().parents[1] / 'scripts' / 'analyze_confirm.py'
sys.argv = [str(script)]
runpy.run_path(str(script), run_name='__main__')
