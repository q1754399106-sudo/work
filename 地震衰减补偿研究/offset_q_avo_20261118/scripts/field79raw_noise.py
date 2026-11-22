"""Post-hoc diagnostic (PLAN_FIELD_I_RAW_ZH.md section 4: the raw records also gave a negative Q -> re-diagnose; not a verdict).
Hypothesis: an ambient-noise floor flattens the far-offset spectra at high frequency (far traces are weaker).
Noise: the raw records are not muted, so traces with 3000 <= |x| <= 4500 m have a pre-arrival window 0.1-0.5 s (water direct
wave at |x| / 1500 >= 2 s, refractions later than 1 s at these offsets). Mean power of that window (same 0.4 s Hann) per
block and component = noise power N(f) of any 0.4 s window (stationary ambient noise, raw records have no gain).
1. SNR P / N at 20, 40, 60 Hz for the near reference and the far classes, windows 1.2-2.2 s.
2. Estimate again with the frozen synthetic estimator's noise rule (oqa.c_off): subtract N from both powers and use only the
   band where both exceed 3 N (contiguous from 12 Hz, >= 10 Hz wide), otherwise as field79raw (pz, a, b).
-> results/field79raw_noise.json"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import field79 as F, field79raw as W  # noqa: E402

NWIN = (0.1, 0.5); NOFF = (3000.0, 4500.0)


def load():
    H, h = W.headers(); sel = (h["dline"] == 5081) & (np.abs(h["cmp"]) < F.CMP_LIM)
    key = lambda i: (h["ffid"][i], h["dpt"][i]); ia = {key(i): i for i in np.where(sel & (h["trid"] == 11))[0]}
    ib = {key(i): i for i in np.where(sel & (h["trid"] == 12))[0]}; pairs = sorted(set(ia) & set(ib))
    ka = np.array([ia[p] for p in pairs]); kb = np.array([ib[p] for p in pairs]); cmp = h["cmp"][ka]; off = h["off"][ka]
    zw = float(np.median(h["wd"][ka])); e = np.linspace(-F.CMP_LIM, F.CMP_LIM, F.NBLK + 1); blocks = []
    for i in range(F.NBLK):
        m = (cmp >= e[i]) & (cmp < e[i + 1]); blocks.append((W.samples(H, ka[m]), W.samples(H, kb[m]), off[m]))
    return blocks, zw


def noise(A, off):
    k = (off >= NOFF[0]) & (off <= NOFF[1]); return F.power(A[k], *NWIN)


def rows_noise(A, B, off, tabs, ptabs, zw, c0sq, comp, NA, NB):
    rows = []; kr = off <= F.REF; xr = float(np.sqrt(np.mean(off[kr] ** 2)))
    comb = lambda PA, PB, cw: PA + PB / (c0sq * cw ** 2) if comp == "pz" else (PA if comp == "a" else PB)
    ib = (F.F >= F.BAND[0]) & (F.F <= F.BAND[1]); snr = []
    for tc, tab, ptab in zip(F.CENT, tabs, ptabs):
        r0 = tab(np.array([xr])); tr = float(r0[0][0]); cwr, taur = W.ghost(ptab(xr)[0], zw)
        PA0, PB0 = F.power(A[kr], tr - F.WL / 2, tr + F.WL / 2), F.power(B[kr], tr - F.WL / 2, tr + F.WL / 2)
        if PA0 is None or PB0 is None:
            continue
        R0 = comb(PA0, PB0, cwr); N0 = comb(NA, NB, cwr)
        for o in F.OFFS:
            k = (off >= o - F.HW) & (off < o + F.HW)
            if k.sum() < 10:
                continue
            xo = float(np.sqrt(np.mean(off[k] ** 2))); tx, la, lb, ls = (float(v[0]) for v in tab(np.array([xo])))
            if not np.isfinite(tx) or tx / tc > F.STRETCH or tx - F.WL / 2 < F.MUTE(o + F.HW) + F.MUTE_PAD or tx + F.WL / 2 > W.T[-1]:
                continue
            PA, PB = F.power(A[k], tx - F.WL / 2, tx + F.WL / 2), F.power(B[k], tx - F.WL / 2, tx + F.WL / 2)
            if PA is None or PB is None:
                continue
            cw, tau = W.ghost(ptab(xo)[0], zw); Rx = comb(PA, PB, cw); Nx = comb(NA, NB, cw)
            snr.append(dict(tc=float(tc), off=o, near=[float(np.interp(f, F.F, R0 / N0)) for f in (20, 40, 60)], far=[float(np.interp(f, F.F, Rx / Nx)) for f in (20, 40, 60)]))
            ok = ib & (R0 > 3 * N0) & (Rx > 3 * Nx); i0 = np.where(F.F >= F.BAND[0])[0][0]; j = i0
            while j < len(F.F) and ok[j]:
                j += 1
            if F.F[j - 1] - F.F[i0] < F.MINBW:
                continue
            s = slice(i0, j); P0 = W.smooth_period(np.maximum(R0 - N0, 1e-12 * R0.max()), taur); Px = W.smooth_period(np.maximum(Rx - Nx, 1e-12 * Rx.max()), tau)
            y = float(-np.polyfit(F.F[s], 0.5 * np.log(Px[s] / P0[s]), 1)[0])
            rows.append(dict(tc=float(tc), off=xo, y=y, fmax=float(F.F[j - 1]), dA=la - float(r0[1][0]), dB=lb - float(r0[2][0]), dS=ls - float(r0[3][0])))
    return rows, snr


def main():
    blocks, zw = load(); R = F.Rays(1.0, zw); tabs = [R.table(tc) for tc in F.CENT]; ptabs = [W.ptable(R, tc) for tc in F.CENT]
    c0sq = W.c0sq_near(blocks, tabs); out = dict(c0sq=c0sq, water_depth=zw)
    for comp in ("pz", "a", "b"):
        rows, snr = [], []
        for A, B, off in blocks:
            r, s = rows_noise(A, B, off, tabs, ptabs, zw, c0sq, comp, noise(A, off), noise(B, off)); rows.append(r); snr += s
        S = F.summarise(rows); far = np.array([q["far"] for q in snr if q["tc"] >= 1.2]); near = np.array([q["near"] for q in snr if q["tc"] >= 1.2])
        S["snr_median_20_40_60"] = dict(near=list(np.median(near, 0)), far_by_offset={str(o): list(np.median([q["far"] for q in snr if q["tc"] >= 1.2 and q["off"] == o], 0)) for o in F.OFFS if any(q["off"] == o and q["tc"] >= 1.2 for q in snr)})
        out[comp] = S
        print(comp, json.dumps({k: S[k] for k in ("invQ_eff", "ci90_eff", "per_block_invQ_eff", "n_rows")}), flush=True)
        print("   SNR near (20/40/60 Hz):", np.round(np.median(near, 0), 1), "| far 1500 / 2500 / 3500 m:",
              [np.round(S["snr_median_20_40_60"]["far_by_offset"].get(str(o), [np.nan] * 3), 1).tolist() for o in (1500, 2500, 3500)], flush=True)
        print("   windows", {k: (v["n"], round(v["inv_q_eff_x1000"], 1), round(v["fmax_median"], 1)) for k, v in S["q_eff_by_window"].items()}, flush=True)
    (F.ROOT / "results" / "field79raw_noise.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
