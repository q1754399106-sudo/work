"""analysis/field_bands.json: shallow / target band edges of lib.wavelets and of the notch-tolerant variant on
the experiment-13 real inlines (whole inlines, xlines 2100-2500 every 2, 1.2-3.0 s). Data-level only: no
compensation is run. The experiment-13 T2 scores of SC-a and IQ70 are attached for a descriptive comparison."""
import sys, json, subprocess
from pathlib import Path
from libc import ROOT, BASE, wavelets, np

E13 = BASE / "field_baselines_volve_20261029"
sys.path.insert(0, str(E13 / "scripts"))
import lib13 as L  # noqa: E402

cm = L.semi.curtain_cm(); lines = json.loads((E13 / "analysis" / "lines.json").read_text(encoding="utf-8"))["lines"]


def e13(ln, tag):
    p = f"地震衰减补偿研究/field_baselines_volve_20261029/results/lines/IL_{ln}_{tag}.json"
    r = subprocess.run(["git", "show", f"HEAD:{p}"], capture_output=True, cwd=str(BASE.parent))
    return json.loads(r.stdout.decode("utf-8")) if r.returncode == 0 else None


rows = {}
for ln in lines:
    D = L.read_inline(ln); w0 = wavelets(D, cm, False); w1 = wavelets(D, cm, True)
    a = e13(ln, "A"); b = e13(ln, "B")
    rows[ln] = dict(fs=w0["fmax_shallow"], fs_robust=w1["fmax_shallow"], fd=w0["fmax_target"], f_peak=w0["f_peak"],
                    live=int((np.abs(D).sum(axis=0) > 0).sum()), e13_A=a, e13_B=b)
    print(ln, {k: v for k, v in rows[ln].items() if not k.startswith("e13")}, flush=True)
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "field_bands.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
