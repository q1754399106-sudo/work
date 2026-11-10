"""Exploratory check of the depth-stationary colour model (colour.py) on the used FD data H1-H6: attenuation
multiplier u (truth 1, colour-free reference 0.96-1.00) from noise-free data and at 30 / 14 dB.

  free   C ~ alpha_w beta(f) |G_u|^2                      colour stationary in time (spectral-ratio assumption)
  depth  C ~ alpha_w beta(f) Phi(f / v_w) |G_u|^2         colour stationary in depth
  truth  C ~ alpha_w T_w(f) |G_u|^2                       window spectra of the elastic target (not available in practice)

usage: pilot_colour.py [n_knots]
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, colour  # noqa: E402
from uqp import np, lib, ROOT, UGRID  # noqa: E402
import fdcheck as fd  # noqa: E402

UG = UGRID[::2]; NK = int(sys.argv[1]) if len(sys.argv) > 1 else 12; out = {}
for name in fd.MODELS:
    D0, T, cm, tw, fp = fd.load_fd(name); x = lib.absorption(cm, tw); rec = {}
    PT = uqp.win_spectra(T, tw)[4]
    for tag, db in (("clean", None), ("30dB", 30.0), ("14dB", 14.0)):
        D = D0 if db is None else D0 + fd.noise(D0, db, 20261122000 + 1000 * (fd.MODELS.index(name) + 1))
        f, cen, C, N, P = uqp.win_spectra(D, tw); Cx = P if db is None else C; wt = uqp.weights(f, P)
        tab = uqp.atten_table(tw, x, x, f, cen, UG); v_w = colour.window_velocity(tw, colour.layered_velocity(cm, tw), cen)
        a = uqp.surface_fit(Cx, wt, tab, ugrid=UG); b = colour.surface_fit3(Cx, wt, tab, f, v_w, ugrid=UG, n_knots=NK)
        c = uqp.surface_fit(Cx, wt, tab, ugrid=UG, g=np.log(PT))
        rec[tag] = dict(free=a["u"], depth=b["u"], truth=c["u"], cost_free=a["cost"], cost_depth=b["cost"], phi=b["phi"].tolist())
        print(name, "%-5s" % tag, "u: free %.3f  depth %.3f  truth-colour %.3f   cost free %.4f depth %.4f" %
              (a["u"], b["u"], c["u"], a["cost"], b["cost"]), " v_w", np.round(v_w[[0, len(v_w) // 2, -1]]).astype(int), flush=True)
    out[name] = rec
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / f"pilot_colour_k{NK}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
