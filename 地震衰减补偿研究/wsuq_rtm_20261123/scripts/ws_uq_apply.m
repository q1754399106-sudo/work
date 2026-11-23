function [R0, UQW1_old, WS_UQ0, WS_UQ1] = ws_uq_apply(Cs, TS, f, sig, fref, iv, s70, fv_band, shot_w)
% Multi-shot Matched-Filter Bayesian UQ-RTM operator (WS-UQ0, WS-UQ1) alongside exact UQW1_old & R0.
% Inputs:
%   Cs      : [nin x ns] complex single/double imaging product conj(Us).*V at frequency f
%   TS      : [nin x ns] real two-way t*/m (already divided by m!) at each interior node and shot
%   f       : scalar frequency (Hz)
%   sig     : log-normal prior std on u (0.447)
%   fref    : reference frequency (Hz); if fv_band = [f1 f2] is non-empty, uses VA-UQ sinc expectation
%   iv      : [nin x 1] stacked 1/snr estimate
%   s70     : stabilization floor (gain_s(70))
%   fv_band : [] for fixed fref, or [f1 f2] for uniform ln(fv) in [f1, f2] (VA-UQ)
%   shot_w  : optional [nin x ns] nonnegative weights (default uniform 1/ns)
%
% Exact boundary identities:
%   1) TS == 0  =>  WS_UQ1 == R0
%   2) TS(:,s) identical across s (and uniform shot_w) => WS_UQ1 == UQW1_old (bit-for-bit to ~1e-6)
%   3) Heterogeneous TS => weights each shot by E[G_s^*] ~ exp(-pi f t_s^*), suppressing dead/edge shots.

[nin, ns] = size(Cs);
if nargin < 8, fv_band = []; end
if nargin < 9 || isempty(shot_w)
    a = ones(nin, ns, 'double') / ns;
else
    a = double(shot_w);
    a = a ./ max(sum(a, 2), realmin);
end

Cs_d = double(Cs);
TS_d = double(TS);
iv_d = double(iv(:));

% Quadrature matching uq_tab.m (41-pt Gauss-Hermite via Golub-Welsch on standard normal)
[z_gh, w_gh] = gh_std_normal(41);
u_gh = exp(sig * z_gh);  % 1 x J

% Dispersion kernel k(f)
use_va = ~isempty(fv_band);
if use_va
    f1 = fv_band(1); f2 = fv_band(2);
    fbar = sqrt(f1 * f2);
    Lam = log(f2 / f1);
    kf = 2 * f * log(max(f, realmin) / fbar);
else
    kf = 2 * f * log(max(f, realmin) / fref);
end
af = -pi * f;

% Build 1D lookup table on xg spanning [0, max(TS)] for fast exact E[G], E|G|^2 (zero 3D memory!)
tmax = max(TS_d(:));
if tmax <= 0
    R0 = real(sum(Cs_d, 2));
    UQW1_old = R0; WS_UQ0 = R0; WS_UQ1 = R0;
    return;
end
ng = 1201;
xg = linspace(0, 1.0001 * tmax, ng)';   % ng x 1
dxg = xg(2) - xg(1);
x_u = xg * u_gh;                         % ng x J
G_tab = exp((af + 1i * kf) * x_u);       % ng x J
if use_va
    arg = f * x_u * Lam;
    sc = ones(size(arg));
    nz = abs(arg) > 1e-12;
    sc(nz) = sin(arg(nz)) ./ arg(nz);
    Gbar_tab = (G_tab .* sc) * w_gh(:);  % E_{u,fv}[G]
else
    Gbar_tab = G_tab * w_gh(:);          % E_u[G]
end
Eg2_tab = (abs(G_tab) .^ 2) * w_gh(:);   % E_u[|G|^2] (same for VA since |e^{i*phase}|=1)

% Linear interpolation index for TS_d (nin x ns)
pos = TS_d / dxg;
i0 = min(max(floor(pos), 0), ng - 2);
ww = pos - i0;
i0 = i0 + 1;

gbar = Gbar_tab(i0) .* (1 - ww) + Gbar_tab(i0 + 1) .* ww;  % nin x ns = E[G_s]
eg2  = Eg2_tab(i0)  .* (1 - ww) + Eg2_tab(i0 + 1)  .* ww;  % nin x ns = E|G_s|^2

% --- R0 ---
R0 = real(sum(Cs_d, 2));

% --- Exact UQW1_old (migrate2.m line 72) ---
A1 = conj(gbar) ./ max(eg2, realmin);
C1 = 1 ./ max(eg2, realmin);
P1 = conj(gbar) ./ max(abs(gbar), realmin);
P1(abs(gbar) == 0) = 1;
H_old = (A1 + s70 * P1 .* C1 + P1 .* iv_d) ./ (1 + s70 * C1 + iv_d);
UQW1_old = real(sum(H_old .* Cs_d, 2));

% --- WS-UQ0 & WS-UQ1 (O(J * nin * ns) scalar loop, zero 3D allocation!) ---
E_gbar = sum(a .* (abs(gbar) .^ 2), 2);     % nin x 1
L_gbar = sum(a .* abs(gbar), 2);            % nin x 1
D_sig = zeros(nin, 1);
D_r0  = zeros(nin, 1);

% Note: for VA-UQ, f_v is global across shots for a given imaging point!
% Under uniform ln(f_v), E_{u,fv}[|sum a_s gbar_s^* G_s(u, fv)|^2] uses exact Gauss-Legendre in fv (5 pts) or center fbar.
% Check: at fixed u, |sum a_s gbar_s^* exp((-pi f + i 2f ln(f/fv)) t_s^* u)|^2 varies very weakly with fv (first zero > 350 Hz),
% using 5-point Gauss-Legendre on ln(fv) in [ln f1, ln f2] is exact to < 1e-7!
if use_va
    [fv_nodes, fv_w] = gl_nodes(5, log(f1), log(f2));
    fv_w = fv_w / sum(fv_w);
    kf_list = 2 * f * (log(max(f, realmin)) - fv_nodes);
else
    kf_list = kf;
    fv_w = 1.0;
end

gbar_conj_a = a .* conj(gbar);              % nin x ns
for j = 1:numel(u_gh)
    uj = u_gh(j);
    wj = w_gh(j);
    amp_j = exp(af * TS_d * uj);            % nin x ns (real!)
    for m = 1:numel(kf_list)
        Gj = amp_j .* exp((1i * kf_list(m)) * (TS_d * uj));
        s_sig = sum(gbar_conj_a .* Gj, 2);  % nin x 1
        s_r0  = sum(a .* Gj, 2);            % nin x 1
        wm = wj * fv_w(m);
        D_sig = D_sig + wm * (abs(s_sig) .^ 2);
        D_r0  = D_r0  + wm * (abs(s_r0)  .^ 2);
    end
end

matched_stack = sum(conj(gbar) .* Cs_d, 2);  % nin x 1
rho_eff = iv_d .* D_r0 + s70;
denom = D_sig + rho_eff .* E_gbar;
beta0 = E_gbar ./ max(denom, realmin);
beta1 = (E_gbar + rho_eff .* (E_gbar ./ max(L_gbar, realmin))) ./ max(denom, realmin);

WS_UQ0 = real(beta0 .* matched_stack);
WS_UQ1 = real(beta1 .* matched_stack);
end

function [z, w] = gh_std_normal(n)
% n-point Gauss-Hermite nodes and weights for standard normal N(0,1)
i = (1:n-1)';
b = sqrt(i / 2);
J = diag(b, 1) + diag(b, -1);
[V, D] = eig(J);
[x, idx] = sort(diag(D));
w = (V(1, idx) .^ 2);
w = w / sum(w);
z = sqrt(2) * x(:)';
end

function [x, w] = gl_nodes(n, a, b)
i = (1:n-1)';
beta = i ./ sqrt(4 * i.^2 - 1);
J = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(J);
[t, idx] = sort(diag(D));
w = 2 * (V(1, idx) .^ 2);
x = 0.5 * (b - a) * t(:)' + 0.5 * (a + b);
w = 0.5 * (b - a) * w(:)';
end
