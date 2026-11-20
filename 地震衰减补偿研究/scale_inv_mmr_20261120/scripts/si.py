"""Scale-invariant (free-gain) Q compensation operators (DERIVATION_ZH.md).

Well-posed form (candidate A in DERIVATION_ZH.md): for output row t and scenario j the loss of the effective operator
alpha_j H, with a free real gain alpha_j >= 0 per row and scenario, and the noise term scaled with it:
    L_j(H, alpha_j) = sum_f w_f ( |alpha_j H_f G_jf - 1|^2 + rho_jf |alpha_j H_f|^2 ),   w_f = |W_out(f)|^2.
Minimising over alpha_j in closed form, with q_j = sum_f w_f Re(H_f G_jf) and P_j = sum_f w_f |H_f|^2 (|G_jf|^2 + rho_jf):
    alpha_j = max(q_j, 0) / P_j,   min_alpha L_j = E0 - max(q_j, 0)^2 / P_j,   E0 = sum_f w_f,
i.e. 1 - (expected zero-lag correlation)^2 with the target R * W_out (ratio-of-expectations, white reflectivity).
This is scale invariant in H, cross-frequency (non-separable) and its single-scenario optimum is the Wiener filter.

Stationarity for a fixed prior pi (any j-mixture) gives the generalised Bayes form
    H_f = sum_j pi_j alpha_j G*_jf / sum_j pi_j alpha_j^2 (|G_jf|^2 + rho_jf),
so the coupling across frequencies enters only through the J scalars alpha_j(t). The solver alternates the two closed
forms (block-coordinate descent on the biconvex J(H, alpha)); it starts at alpha = 1, i.e. at the ordinary Bayes form
(UQ-W for the lognormal prior), and the objective never increases.

Cells are those of qunknown_comp mmr.cells (normalised by the least-attenuated scenario per (t, f); the substitution
is per frequency and alpha is a scalar, so q, P and every loss are unchanged).

Rules
    SIU  Bayes, pi = UQ-W lognormal Gauss-Hermite weights (sigma 0.447, 41 nodes) - scale-invariant UQ-W
    SIR  minimax regret, one pi for the section over the 17-point scalar family u in [0.5, 2] (as RGJ):
         max_pi mean_t min_{H, alpha} sum_j pi_j R_jt,  R_jt = (min_alpha L_j - L*_j) / E0  (L*_j: Wiener optimum)
         exponentiated gradient on pi (Danskin gradient at the inner local minimum), warm-started alpha;
         compute caps fixed after the smoke run (before any score): pi from every 4th row, at most 150 outer steps,
         at most 20 warm-started alternations per step; the final operator alternates to convergence (<= 100).
"""
import sys
from pathlib import Path

MAIN = Path(r"D:\ChatGPT\研究\地震衰减补偿研究")          # data caches (git-ignored) live only in the main checkout
sys.path.insert(0, str(MAIN / "qunknown_comp_20261109" / "scripts"))
from mmr import np, semi, DT, NFFT, SIG, G_IQ, GH, U1, setup, cells, apply_rows, snr_rows  # noqa: E402

BLK = 64


def wspec(w_out, nfft=NFFT):
    """w_f = |W_out(f)|^2 on the operator grid, normalised to sum 1 (scale of w is irrelevant)."""
    W = np.abs(np.fft.rfft(np.asarray(w_out, float), nfft)) ** 2; return W / W.sum()


class Rows:
    """Pre-computed cells for a set of rows and scenarios; all arrays (n, nf, J)."""
    def __init__(self, xs, f, lf, s, inv, w):
        self.Gc, self.den, ls, self.sc = cells(xs, f, lf, s, inv)
        self.G = np.conj(self.Gc); self.w = w; self.E0 = w.sum()
        self.Lstar = np.matmul(w[None, None, :], ls)[:, 0, :]                    # Wiener optimum per row and scenario

    def h(self, pi, al):
        a = pi[None, :] * al
        return np.matmul(self.Gc, a[..., None].astype(complex))[..., 0] / np.matmul(self.den, (a * al)[..., None])[..., 0]

    def qp(self, H):
        q = np.matmul((self.w * H)[:, None, :], self.G)[:, 0, :].real
        P = np.matmul((self.w * np.abs(H) ** 2)[:, None, :], self.den)[:, 0, :]; return q, P

    def solve(self, pi, al=None, its=100, tol=1e-8):
        """Alternating minimisation of sum_j pi_j L_j(H, alpha_j) for every row; returns H (normalised), alpha, info."""
        n, J = self.Gc.shape[0], self.Gc.shape[-1]; al = np.ones((n, J)) if al is None else al.copy(); obj_prev = None
        for it in range(its):
            H = self.h(pi, al); q, P = self.qp(H); al = np.maximum(q, 0) / P
            obj = (self.E0 - np.maximum(q, 0) ** 2 / P) @ pi                 # (n,) objective after the alpha step
            if obj_prev is not None and np.max(obj_prev - obj) < tol * max(float(obj.max()), 1e-12):
                break
            obj_prev = obj
        k = (al * pi[None, :]).sum(1, keepdims=True); k = np.where(k > 0, k, 1.0)     # level: prior-mean gain = 1
        return H * k, al / k, dict(its=it + 1)

    def regret(self, H):
        q, P = self.qp(H); return (self.E0 - np.maximum(q, 0) ** 2 / P - self.Lstar) / self.E0, q / P   # (n, J)


def solve_sir(R, its=150, tol=1e-5, inner=20):
    """max_pi mean_t min_{H, alpha} sum_j pi_j R_jt; exponentiated gradient with backtracking (as mmr.solve_joint)."""
    J = R.Gc.shape[-1]; pi = np.full(J, 1.0 / J); eta = 10.0; n_in = 0

    def ev(p, al):
        H, al2, inf = R.solve(p, al, its=inner); Rg = R.regret(H)[0].mean(0); return Rg, float(p @ Rg), al2, inf["its"]
    Rg, phi, al, k = ev(pi, None); n_in += k
    for it in range(its):
        q = pi * np.exp(eta * (Rg - Rg.max())); q /= q.sum(); Rq, pq, alq, k = ev(q, al); n_in += k
        if pq >= phi - 1e-15:
            pi, Rg, phi, al = q, Rq, pq, alq; eta = min(eta * 1.5, 1e6)
        else:
            eta *= 0.3
        if Rg.max() - phi < tol * max(abs(phi), 1e-12) or eta < 1e-8:
            break
    return pi, dict(gap=float(Rg.max() - phi), value=float(phi), its=it + 1, inner=n_in, regret=Rg.tolist())


def scenarios(x_hat, dx_hat, rule):
    if rule == "SIU":
        e, w = GH; return x_hat[:, None] + dx_hat[:, None] * (np.exp(SIG * e) - 1)[None, :], w / w.sum()
    return x_hat[:, None] + dx_hat[:, None] * (U1 - 1)[None, :], None


def design(D, cm, Q, c, tw, t_anchor, rule, inv, w, g_db=G_IQ, sub=4):
    """Output of one rule and info. SIR: pi from every `sub`-th row; operator on all rows (in blocks)."""
    f, lf, x_hat, dx_hat = setup(cm, Q, c, tw, t_anchor); s = semi.gain_sigma2(g_db); xs, pi = scenarios(x_hat, dx_hat, rule)
    info = {}
    if rule == "SIR":
        pi, info = solve_sir(Rows(xs[::sub], f, lf, s, inv[::sub], w)); info["pi"] = pi.tolist()
    Hs, al_all, its = [], [], []
    for a in range(0, len(xs), BLK):
        b = slice(a, a + BLK); R = Rows(xs[b], f, lf, s, inv[b], w); H, al, inf = R.solve(pi)
        Hs.append(H * R.sc); al_all.append(al); its.append(inf["its"])
    al = np.concatenate(al_all); info.update(alt_its_max=int(max(its)), alpha_q=np.quantile(al, [0.05, 0.5, 0.95], axis=0).tolist())
    return apply_rows(np.concatenate(Hs), D), info


def bayes_rows(x_hat, dx_hat, rule):
    """For checks: scenario profiles and prior of a rule."""
    return scenarios(x_hat, dx_hat, rule)
