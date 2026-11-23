function posthoc_g2_dip_check()
% Post-hoc mechanism verification ONLY (does NOT alter pre-registered spg_verdict.json!):
% Tests whether replacing flat-mirror x_mir = 2x - x_s with Dip-Aware Specular Mirror
% x_mir(phi) = x + z * tan(theta_s + 2*phi) (with local structural dip phi estimated strictly
% blindly from the uncompensated R0 image via 2D gradient structure tensor) eliminates the
% -1.37 dB SNR & -0.025 E1_HF loss observed on steep 18-deg dipping model G2!
maxNumCompThreads(2);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out_dir = fullfile(root, 'results');
addpath('D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\scripts');

t0_all = tic;
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx;
P.x = ((1:P.nx)  - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)};
P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:6:152);
P.ixr = P.npml + 1 + (0:2:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2});
[Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

% Exact G2 model (seed = 20261132)
seed = 20261132; rng(seed, 'twister');
Zg = Zfull;
hw = (Xfull > 760.0 + (Zfull - 200.0) / tand(68.0)) & (Zfull > 200.0);
Zg(hw) = Zg(hw) - 35.0;
z1 = 240.0 + 0.08 * (Xfull - 800.0);
z2 = 650.0 + 0.32 * (Xfull - 800.0);
z3 = 790.0 - 0.25 * (Xfull - 800.0);
z4 = 900.0 + 0.10 * (Xfull - 800.0);
v_true = 1850.0 * ones(P.nz, P.nx);
v_true(Zg >= z1) = 2150.0; v_true(Zg >= z2) = 2600.0;
v_true(Zg >= z3) = 2950.0; v_true(Zg >= z4) = 3320.0;
v_mig = smooth2d(v_true, 2.5);
invQ = (1/85)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
lens = exp(-(((Xfull - 560)/210).^2 + ((Zfull - 480)/90).^2).^2);
invQ = invQ + (1/18 - 1/85)*lens;
Q_true = 1./invQ;

% Eikonal fields (keep F.ts for dip-aware two_way!)
s_mig = 1./v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k) = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns);
px_s = zeros(nin, ns); pz_s = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    Tk = double(F.T(P.in{1}, P.in{2}, k));
    [gx, gz] = gradient(Tk, P.dx); gn = max(sqrt(gx.^2 + gz.^2), 1e-12);
    px_s(:, k) = reshape(gx ./ gn, nin, 1); pz_s(:, k) = reshape(gz ./ gn, nin, 1);
end

freqs = 4.0:1.0:37.0; nf = numel(freqs);
fref = 25.0; fpeak = 22.0; sig = 0.447; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
flip_sgn = 2 * (rand(1, ns) > 0.5) - 1;

% Calibrate noise identically to spg_rtm_bench
W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = Uv22(irec,:) - Dbg22; clear Uv22 Lv Uv Pv Qv Rv
d_rms = sqrt(mean(abs(Dcal(:)).^2));
sig_noisy = d_rms * 10^(-3/20);

% Pass 1 (cache single-freq Cs_sig, Cs_nu, I_ref in RAM: 34*16000*25*8B = 108 MB -> super fast!)
Cs_sig_all = zeros(nin, ns, nf, 'single');
Cs_nu_all  = zeros(nin, ns, nf, 'single');
I_ref = zeros(nin, 1);
I_R0_noisy = zeros(nin, 1);
for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
    [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
    U_bg = Qm*(Um\(Lm\(Pm*(Rm\b_src)))); D_bg = U_bg(irec,:);
    Us_c = conj(U_bg(lin,:)) * w^2; clear U_bg
    [Le,Ue,Pe,Qe,Re] = lu(helm_op(v_true, f, P));
    U_el = Qe*(Ue\(Le\(Pe*(Re\b_src)))); D_el = U_el(irec,:) - D_bg; clear U_el Le Ue Pe Qe Re
    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = U_vi(irec,:) - D_bg; clear U_vi Lv Uv Pv Qv Rv b_src
    N_unit = (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns));
    bp = @(D) extract_bp(D, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_el = bp(D_el); V_vi = bp(D_vi); V_nu = bp(N_unit); clear Lm Um Pm Qm Rm
    I_ref = I_ref + real(sum(Us_c .* V_el, 2)); clear V_el
    cs = Us_c .* V_vi; cn = Us_c .* V_nu;
    Cs_sig_all(:,:,jf) = single(cs);
    Cs_nu_all(:,:,jf)  = single(cn);
    I_R0_noisy = I_R0_noisy + real(sum(cs + sig_noisy*cn, 2));
end

% Estimate local structural dip phi(z,x) BLINDLY from noisy R0 via 2D gradient structure tensor!
R0_2d = prep_img(I_R0_noisy, P.nzi, P.nxi, z_in);
[gx, gz] = gradient(R0_2d, P.dx);
Jxx = smooth2d(gx.^2, 4.0); Jzz = smooth2d(gz.^2, 4.0); Jxz = smooth2d(gx.*gz, 4.0);
% Dominant eigenvector (normal to reflector) -> angle from vertical!
% Normal (nx, nz) = (sin(phi), cos(phi)) where dip slope dz/dx = -nx/nz -> phi = 0.5*atan2(2*Jxz, Jzz - Jxx)
% Wait: for z = z0 + s*x, reflector normal is (-s, 1)/sqrt(1+s^2), so Jxz < 0 when s > 0 -> dip_rad = -0.5*atan2(2*Jxz, Jzz - Jxx)!
% Verify: gx/gz = -s -> 2*Jxz/(Jzz-Jxx) = -2s/(1-s^2) = -tan(2*phi) => phi = -0.5*atan2(2*Jxz, Jzz - Jxx)!
phi_est = -0.5 * atan2(2*Jxz, Jzz - Jxx);
% Clip estimated dip to [-30, +30] deg for stability
phi_est = max(min(phi_est, deg2rad(30)), -deg2rad(30));
fprintf('Blind R0 structure-tensor dip check on G2: r2(true +17.7 deg) med=%.1f deg | r3(true -14.0 deg) med=%.1f deg\n', ...
    rad2deg(median(phi_est(z_in>=600&z_in<=700, x_in>=600&x_in<=1000), 'all')), ...
    rad2deg(median(phi_est(z_in>=760&z_in<=840, x_in>=600&x_in<=1000), 'all')));

% Build Dip-Aware Specular Mask w_dsa AND Dip-Aware two-way TS_dip using x_mir(phi) = x + z*tan(theta_s + 2*phi)!
xr_lo = min(xr); xr_hi = max(xr);
w_sa_flat = zeros(nin, ns);
w_dsa     = zeros(nin, ns);
TS_dip    = zeros(nin, ns);
phi_v = phi_est(:);
Fts_in = zeros(nin, ns);
for k = 1:ns
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
for k = 1:ns
    th_s = atan2(Xin(:) - xs(k), max(Zin(:) - 20.0, 10.0)); % incidence angle from vertical
    xmir_flat = 2*Xin(:) - xs(k);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir_dip = Xin(:) + max(Zin(:) - 20.0, 10.0) .* tan(th_r);

    t_graz = cos_taper(80.0 - abs(rad2deg(th_s)), 0.0, 10.0);
    t_graz_r = cos_taper(80.0 - abs(rad2deg(th_r)), 0.0, 10.0);
    w_sa_flat(:, k) = max(cos_taper(min(xmir_flat - xr_lo, xr_hi - xmir_flat), -80, 120) .* t_graz, 1e-4);
    w_dsa(:, k)     = max(cos_taper(min(xmir_dip  - xr_lo, xr_hi - xmir_dip),  -80, 120) .* t_graz .* t_graz_r, 1e-4);

    % Interpolate receiver-leg t* at xmir_dip!
    q = (xmir_dip - xs(1)) / (xs(2) - xs(1)) + 1;
    q = min(max(q, 1), ns); k0 = floor(q); k0(k0 == ns) = ns - 1; ww = q - k0;
    idx0 = sub2ind([nin ns], (1:nin)', k0); idx1 = sub2ind([nin ns], (1:nin)', k0+1);
    TS_dip(:, k) = Fts_in(:, k) + (1 - ww).*Fts_in(idx0) + ww.*Fts_in(idx1);
end

% Compare RAW_UQW1, SA_flat_UQW1, DSA_UQW1 (dip-aware mask + flat t*), DSA2_UQW1 (dip-aware mask + dip-aware t*) on noisy G2!
xmax = 1.05 * max(max(TS_flat(:)), max(TS_dip(:))) / min(ms);
xg = linspace(0, xmax, 1201); dxg = xg(2);
cfg = {
    'RAW_UQW1',     ones(nin, ns), TS_flat;
    'SA_flat_UQW1', w_sa_flat,     TS_flat;
    'DSA_UQW1',     w_dsa,         TS_flat;
    'DSA2_UQW1',    w_dsa,         TS_dip;
    'RAW_R0',       ones(nin, ns), TS_flat;
    'DSA_R0',       w_dsa,         TS_dip;
};
ncfg = size(cfg, 1);
IS = cell(ncfg, 1); IN = cell(ncfg, 1);
for ic = 1:ncfg
    nc = 1 + (nm-1)*~endsWith(cfg{ic,1}, 'R0');
    IS{ic} = zeros(nin, nc); IN{ic} = zeros(nin, nc);
end

for jf = 1:nf
    f = freqs(jf);
    Cs_s = double(Cs_sig_all(:,:,jf));
    Cs_n = sig_noisy * double(Cs_nu_all(:,:,jf));
    Cs_t = Cs_s + Cs_n;
    T1 = uq_tab(f, xg, sig, fref);
    for ic = 1:ncfg
        Wm = cfg{ic,2}; TSm = cfg{ic,3};
        if endsWith(cfg{ic,1}, 'R0')
            IS{ic} = IS{ic} + real(sum(Wm .* Cs_s, 2));
            IN{ic} = IN{ic} + real(sum(Wm .* Cs_n, 2));
            continue;
        end
        stk_t = abs(sum(Wm .* Cs_t, 2)).^2;
        stk_f = abs(sum(Wm .* Cs_t .* flip_sgn, 2)).^2;
        pt = smooth2d(reshape(stk_t, P.nzi, P.nxi), 5.0);
        pf = smooth2d(reshape(stk_f, P.nzi, P.nxi), 5.0);
        iv = min(max(pf(:) ./ max(pt(:) - pf(:), 0.05*pf(:)), 1e-3), 1e3);
        for im = 1:nm
            u = (TSm / ms(im)) / dxg; i0 = min(floor(u), numel(xg)-2); ww = u - i0; i0 = i0 + 1;
            it = @(tb) tb(i0).*(1-ww) + tb(i0+1).*ww;
            A1 = it(T1.A); C1 = it(T1.C); P1 = it(T1.P);
            Huq = (A1 + s70*P1.*C1 + P1.*iv) ./ (1 + s70*C1 + iv);
            IS{ic}(:, im) = IS{ic}(:, im) + real(sum((Wm.*Huq).*Cs_s, 2));
            IN{ic}(:, im) = IN{ic}(:, im) + real(sum((Wm.*Huq).*Cs_n, 2));
        end
    end
end

ref_2d = prep_img(I_ref, P.nzi, P.nxi, z_in);
ref_hf = highpass_kz(ref_2d, P.dx, 10.0);
iz_w = find(z_in >= 560 & z_in <= 950);
ix_w = find(x_in >= 180 & x_in <= 1420);
fprintf('\n=== Post-hoc G2 (18-deg Steep Dip) Dip-Aware Specular Mask Check [noisy W_deep] ===\n');
out = struct('elapsed_s', toc(t0_all));
for ic = 1:ncfg
    nc = size(IS{ic}, 2); e1 = zeros(1,nc); e1h = zeros(1,nc); snr = zeros(1,nc);
    for c = 1:nc
        S2 = prep_img(IS{ic}(:,c), P.nzi, P.nxi, z_in);
        N2 = prep_img(IN{ic}(:,c), P.nzi, P.nxi, z_in);
        T2 = S2 + N2; Th = highpass_kz(T2, P.dx, 10.0);
        e1(c) = score_e1(T2, ref_2d, iz_w, ix_w);
        e1h(c) = score_e1(Th, ref_hf, iz_w, ix_w);
        snr(c) = 10*log10(sum(S2(iz_w,ix_w).^2,'all') / max(sum(N2(iz_w,ix_w).^2,'all'), 1e-30));
    end
    fprintf('  %-13s | E1=%.4f (m0.5=%.4f,m1=%.4f) | HF=%.4f | SNR=%.2f dB\n', ...
        cfg{ic,1}, mean(e1), e1(1), e1(min(3,nc)), mean(e1h), mean(snr));
    out.(cfg{ic,1}) = struct('E1_mean', mean(e1), 'E1_HF_mean', mean(e1h), 'SNR_dB_mean', mean(snr));
end
fn = fullfile(out_dir, 'posthoc_g2_dip_check.json');
fid = fopen(fn, 'w'); fwrite(fid, jsonencode(out, 'PrettyPrint', true), 'char'); fclose(fid);
fprintf('Saved %s (%.1fs)\n', fn, toc(t0_all));
end

function t = cos_taper(x, x0, x1), u = min(max((x - x0)/max(x1 - x0, 1e-12), 0), 1); t = 0.5 - 0.5*cos(pi*u); end
function V = extract_bp(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns); b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2; V = Vf(lin, :);
end
function B = smooth2d(A, sig)
r = max(1, ceil(3*sig)); g = exp(-0.5*((-r:r)/sig).^2); g = g/sum(g);
Ap = [repmat(A(1,:), r, 1); A; repmat(A(end,:), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:,1), 1, r), B1, repmat(B1(:,end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end
function y = prep_img(v, nzi, nxi, z_in), y = reshape(double(v), nzi, nxi); y = y - smooth2d(y, 5.0); y(z_in < 160, :) = 0; end
function yh = highpass_kz(y, dx, kz_cut)
n = 2^nextpow2(2*size(y,1)); Y = fft(y, n, 1);
kz = abs([0:n/2, -n/2+1:-1]')/(n*dx)*1000; g = 0.5 - 0.5*cos(pi*min(max((kz-(kz_cut-2))/4,0),1));
yh = real(ifft(Y.*g, n, 1)); yh = yh(1:size(y,1), :);
end
function s = score_e1(a, r, iz, ix)
rr = r(iz,ix); c0 = cc(a(iz,ix), rr); cl = c0; nz = size(a,1);
for L = [-2 -1 1 2], cl = max(cl, cc(a(min(max(iz+L,1),nz), ix), rr)); end
s = 0.5*(c0 + cl);
end
function c = cc(a, b), a = a(:)-mean(a(:)); b = b(:)-mean(b(:)); c = (a'*b)/max(sqrt((a'*a)*(b'*b)),1e-30); end
