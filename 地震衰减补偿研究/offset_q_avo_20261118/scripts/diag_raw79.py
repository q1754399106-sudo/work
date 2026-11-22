import sys, numpy as np
sys.path.insert(0, "D:/ChatGPT/研究/.claude/worktrees/study-research-docs-plan-c4d4c6/地震衰减补偿研究/offset_q_avo_20261118/scripts")
import field79 as F, field79raw as W
from multiprocessing import Pool
def job(args):
    g, rsb, sw, nz, ib = args; W.GHOST, W.R_SB, W.SWELL, W.NOISE = g, rsb, sw, nz; return W.synth_block(ib)
if __name__ == "__main__":
    H, h = F.read_headers(); sel = np.abs(h["cmp"]) < F.CMP_LIM; zw = float(np.median(h["wd"][sel])); R = F.Rays(1.0, zw)
    tabs = [R.table(tc) for tc in F.CENT]; ptabs = [W.ptable(R, tc) for tc in F.CENT]
    cases = [(True, 0.3, True, False), (True, 0.3, False, True), (True, 0.3, True, True)]
    with Pool(4) as pool:
        for c in cases:
            blocks = pool.map(job, [c + (ib,) for ib in (0, 1)])
            c0 = W.c0sq_near(blocks, tabs); out = []
            for comp in ("pz", "a", "b"):
                rows = [r for A, B, o in blocks for r in W.rows_pz(A, B, o, tabs, ptabs, zw, c0, comp)]
                x = np.pi * np.array([r["dS"] for r in rows]); yt = np.pi * np.array([r["dA"] / F.QS[0] + r["dB"] / F.QS[1] + (r["dS"] - r["dA"] - r["dB"]) / F.QS[2] for r in rows])
                out.append("%s %.2f (truth %.2f, n %d)" % (comp, F.fit1(rows), 1000 * (x @ yt) / (x @ x), len(rows)))
            # without period smoothing, pz
            print("ghost %s R %.1f swell %s noise %s  c0sq %.3f | %s" % (c[0], c[1], c[2], c[3], c0, " | ".join(out)), flush=True)
