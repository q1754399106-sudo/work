function run_unified_synth_midSNR()
% Calibration-corrected supplementary synthetic suite (U1b-U4b) at realistic +2 to +5 dB image SNR.
% Preserves pre-registered b6018cf run (synth_U1..U5.json, which exposed the -14 dB blind-dip floor)
% completely untouched, and evaluates the exact same U1-U4 models/seeds at realistic l30 image SNR
% (snr_db = -2.0 dB for U1b-U3b, +2.0 dB for U4b -> ~+2 to +5 dB W_deep image SNR).

maxNumCompThreads(8);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res_dir = fullfile(root, 'results');
wt_root = fileparts(root);
addpath(fullfile(wt_root, 'uq_rtm_20261117', 'scripts'));
core = uq_rtm_unified_core();
cfg  = core.default_cfg('synth_10m');

for mcell = {'U1', 'U2', 'U3', 'U4'}
    mname = mcell{1};
    fn = fullfile(res_dir, ['synth_' mname 'b.json']);
    if isfile(fn)
        fprintf('Skipping %sb (exists)\n', mname); continue;
    end
    if strcmp(mname, 'U4')
        res = run_u4b(core, cfg);
    else
        res = run_ub(mname, core, cfg);
    end
    fid = fopen(fn, 'w'); fwrite(fid, jsonencode(res, 'PrettyPrint', true), 'char'); fclose(fid);
    fprintf('Saved %s (%.1fs)\n\n', fn, res.elapsed_s);
end
end

function res = run_ub(mname, core, cfg)
t0_all = tic;
P.dx = 10.0; P.npml = 18; P.nzi = 100; P.nxi = 160;
P.nz = P.nzi + 2*P.npml; P.nx = P.nxi + 2*P.npml;
P.z = ((1:P.nz)' - P.npml - 1) * P.dx; P.x = ((1:P.nx) - P.npml - 1) * P.dx;
P.in = {P.npml + (1:P.nzi), P.npml + (1:P.nxi)}; P.izs = P.npml + 3;
P.ixs = P.npml + 1 + (8:3:152); P.ixr = P.npml + 1 + (0:1:P.nxi-1);
ns = numel(P.ixs); nr = numel(P.ixr); N = P.nz * P.nx;
[Zi, Xi] = ndgrid(P.in{1}, P.in{2}); lin = sub2ind([P.nz P.nx], Zi(:), Xi(:)); nin = numel(lin);
isrc = sub2ind([P.nz P.nx], P.izs*ones(1,ns), P.ixs);
irec = sub2ind([P.nz P.nx], P.izs*ones(1,nr), P.ixr);
[Zfull, Xfull] = ndgrid(P.z, P.x);
z_in = P.z(P.in{1}); x_in = P.x(P.in{2}); [Zin, Xin] = ndgrid(z_in, x_in);
xs = P.x(P.ixs); xr = P.x(P.ixr);

rec_mask = true(nr, ns); corr_noise_sig = 5.0; snr_db = -2.0;
switch mname
    case 'U1'
        seed = 20261301;
        z1 = 230.0 - 15.0 * exp(-((Xfull - 760.0)/330.0).^2);
        z2 = 625.0 - 38.0 * exp(-((Xfull - 800.0)/300.0).^2) + 0.08*(Xfull - 800.0);
        z3 = 760.0 - 28.0 * exp(-((Xfull - 800.0)/330.0).^2) + 0.06*(Xfull - 800.0);
        z4 = 890.0 - 14.0 * exp(-((Xfull - 800.0)/370.0).^2);
        v_true = 1850.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2160; v_true(Zfull>=z2)=2580; v_true(Zfull>=z3)=2940; v_true(Zfull>=z4)=3320;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        invQ = invQ + (1/14 - 1/90)*exp(-(((Xfull - 790)/210).^2 + ((Zfull - 415)/140).^2).^1.5);
        Q_true = 1./invQ; W_targ = {[580 1000], [560 950]}; W_flank = {[200 480], [560 950]};
    case 'U2'
        seed = 20261302; Zg = Zfull;
        hw = (Xfull < 820.0 - (Zfull - 200.0) / tand(63.0)) & (Zfull > 200.0); Zg(hw) = Zg(hw) + 32.0;
        z1 = 235.0 - 0.05*(Xfull-800); z2 = 655.0 + 0.40*(Xfull-800); z3 = 785.0 - 0.29*(Xfull-800); z4 = 900.0 - 0.10*(Xfull-800);
        v_true = 1860.0 * ones(P.nz, P.nx); v_true(Zg>=z1)=2170; v_true(Zg>=z2)=2630; v_true(Zg>=z3)=2990; v_true(Zg>=z4)=3360;
        invQ = (1/85)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        invQ = invQ + (1/15 - 1/85)*exp(-(((Xfull - 960)/250).^2 + ((Zfull - 365)/100).^2).^1.5);
        Q_true = 1./invQ; W_targ = {[780 1200], [560 950]}; W_flank = {[220 560], [560 950]};
    case 'U3'
        seed = 20261303;
        for s = 1:ns, rec_mask(:, s) = (xr >= xs(s) - 240.0) & (xr <= xs(s) + 860.0); end
        z1 = 240.0 + 0.04*(Xfull-800); z2 = 640.0 + 0.20*(Xfull-800); z3 = 770.0 + 0.17*(Xfull-800); z4 = 890.0 - 0.04*(Xfull-800);
        v_true = 1840.0 * ones(P.nz, P.nx); v_true(Zfull>=z1)=2140; v_true(Zfull>=z2)=2560; v_true(Zfull>=z3)=2910; v_true(Zfull>=z4)=3280;
        invQ = (1/90)*ones(P.nz, P.nx); invQ(Zfull < 80) = 1/1000;
        invQ = invQ + (1/14 - 1/90)*exp(-(((Xfull - 780)/195).^2 + ((Zfull - 415)/145).^2).^1.5);
        Q_true = 1./invQ; W_targ = {[580 980], [560 950]}; W_flank = {[1040 1380], [560 950]};
end

rng(seed, 'twister');
v_mig = core.smooth2d(v_true, 2.5);
s_mig = 1.0 ./ v_mig; eps_q = 0.5;
F.T = zeros(P.nz, P.nx, ns, 'single'); F.ts = F.T;
for k = 1:ns
    F.T(:,:,k)  = eikonal_fs(s_mig, P.dx, P.izs, P.ixs(k));
    F.ts(:,:,k) = (eikonal_fs(s_mig.*(1+eps_q./Q_true), P.dx, P.izs, P.ixs(k)) - ...
                   eikonal_fs(s_mig.*(1-eps_q./Q_true), P.dx, P.izs, P.ixs(k))) / (2*eps_q);
end
TS_flat = zeros(nin, ns); Fts_in = zeros(nin, ns);
for k = 1:ns
    [~, tt] = two_way(F, k, P); TS_flat(:, k) = double(tt(lin));
    tk = F.ts(:,:,k); Fts_in(:, k) = double(tk(lin));
end
clear F
xr_lo = arrayfun(@(k) min(xr(rec_mask(:,k))), 1:ns);
xr_hi = arrayfun(@(k) max(xr(rec_mask(:,k))), 1:ns);
w_deep_mask = (Zin(:) >= 560 & Zin(:) <= 950) & (Xin(:) >= 180 & Xin(:) <= 1420);

freqs = 4.0:1.0:37.0; nf = numel(freqs); fref_phys = 25.0; fpeak = 22.0; sig = cfg.sig; s70 = gain_s(70);
ms = [0.5, 0.7, 1.0, 1.4, 2.0]; nm = numel(ms);
rng(seed + 1, 'twister');   EPS_B1   = 2*(rand(nr, ns) > 0.5) - 1;
rng(seed + 300, 'twister'); sgn_B300 = 2*(rand(1, ns)  > 0.5) - 1;

W22 = (22/fpeak)^2 * exp(-(22/fpeak)^2);
b22 = zeros(N, ns); b22(sub2ind([N ns], isrc, 1:ns)) = -W22 / P.dx^2;
[Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, 22, P)); Ub = Qm*(Um\(Lm\(Pm*(Rm\b22)))); Dbg22 = Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc22 = v_true ./ (1 - (1./(pi*Q_true))*log(22/fref_phys) - 1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc22, 22, P)); Uv22 = Qv*(Uv\(Lv\(Pv*(Rv\b22)))); Dcal = (Uv22(irec,:) - Dbg22) .* rec_mask; clear Uv22 Lv Uv Pv Qv Rv
sig_n = sqrt(mean(abs(Dcal(rec_mask)).^2)) * 10^(-snr_db/20);
rc = ceil(3 * corr_noise_sig); g_nc = exp(-0.5 * ((-rc:rc)' / corr_noise_sig).^2); g_nc = g_nc / sqrt(sum(g_nc.^2));

Cs_sig_all = zeros(nin, ns, nf, 'single'); Cs_nu_all = Cs_sig_all; Cs_b1_all = Cs_sig_all;
I_ref = zeros(nin, 1); I_R0_tot = zeros(nin, 1);
for jf = 1:nf
    f = freqs(jf); w = 2*pi*f;
    W_src = (f/fpeak)^2 * exp(-(f/fpeak)^2);
    b_src = zeros(N, ns); b_src(sub2ind([N ns], isrc, 1:ns)) = -W_src / P.dx^2;
    [Lm,Um,Pm,Qm,Rm] = lu(helm_op(v_mig, f, P));
    U_bg = Qm*(Um\(Lm\(Pm*(Rm\b_src)))); D_bg = U_bg(irec,:); Us_c = conj(U_bg(lin,:)) * w^2; clear U_bg
    [Le,Ue,Pe,Qe,Re] = lu(helm_op(v_true, f, P));
    U_el = Qe*(Ue\(Le\(Pe*(Re\b_src)))); D_el = (U_el(irec,:) - D_bg) .* rec_mask; clear U_el Le Ue Pe Qe Re
    vc = v_true ./ (1 - (1./(pi*Q_true))*log(max(f,1e-6)/fref_phys) - 1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv] = lu(helm_op(vc, f, P));
    U_vi = Qv*(Uv\(Lv\(Pv*(Rv\b_src)))); D_vi = (U_vi(irec,:) - D_bg) .* rec_mask; clear U_vi Lv Uv Pv Qv Rv b_src
    N_f = conv2(g_nc, 1, sig_n * (1/sqrt(2)) * (randn(nr, ns) + 1i*randn(nr, ns)), 'same') .* rec_mask;
    V_el = extract_bp(D_el, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_vi = extract_bp(D_vi, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_nu = extract_bp(N_f,  Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    V_b1 = extract_bp((D_vi + N_f) .* EPS_B1, Lm, Um, Pm, Qm, Rm, irec, lin, N, ns, P.dx);
    clear Lm Um Pm Qm Rm
    I_ref = I_ref + real(sum(Us_c .* V_el, 2)); clear V_el
    cs = Us_c .* V_vi; cn = Us_c .* V_nu; cb = Us_c .* V_b1;
    Cs_sig_all(:,:,jf) = single(cs); Cs_nu_all(:,:,jf) = single(cn); Cs_b1_all(:,:,jf) = single(cb);
    I_R0_tot = I_R0_tot + real(sum(cs + cn, 2));
end

R0_2d = core.prep_image(I_R0_tot, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
[phi_eff, c_st] = core.blind_dip(R0_2d, P.dx, cfg.dip_sigma_m, cfg.dip_cap_deg);
[w_dsa, TS_dsa] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, phi_eff(:), Fts_in, cfg);
[w_sa,  ~     ] = core.build_dsa_in(Xin(:), Zin(:), 20.0, xs, xr_lo, xr_hi, zeros(nin,1), Fts_in, cfg);
dip_stats = struct('med_abs_phi_deg', median(abs(rad2deg(phi_eff(w_deep_mask)))), ...
                   'p90_abs_phi_deg', prctile(abs(rad2deg(phi_eff(w_deep_mask))), 90), ...
                   'med_coh', median(c_st(w_deep_mask)));

I_raw_s = zeros(nin, nf, 'single'); I_raw_n = I_raw_s; I_raw_t = I_raw_s; I_raw_f1 = I_raw_s; I_raw_f3 = I_raw_s;
I_dsa_s = zeros(nin, nf, 'single'); I_dsa_n = I_dsa_s; I_dsa_t = I_dsa_s; I_dsa_f1 = I_dsa_s; I_dsa_f3 = I_dsa_s;
I_sa_t = zeros(nin, nf, 'single'); I_sa_f3 = I_sa_t;
for jf = 1:nf
    cs = double(Cs_sig_all(:,:,jf)); cn = double(Cs_nu_all(:,:,jf)); ct = cs + cn; cb = double(Cs_b1_all(:,:,jf));
    I_raw_s(:,jf)=single(real(sum(cs,2))); I_raw_n(:,jf)=single(real(sum(cn,2)));
    I_raw_t(:,jf)=single(real(sum(ct,2))); I_raw_f1(:,jf)=single(real(sum(cb,2))); I_raw_f3(:,jf)=single(real(sum(ct.*sgn_B300,2)));
    I_dsa_s(:,jf)=single(real(sum(w_dsa.*cs,2))); I_dsa_n(:,jf)=single(real(sum(w_dsa.*cn,2)));
    I_dsa_t(:,jf)=single(real(sum(w_dsa.*ct,2))); I_dsa_f1(:,jf)=single(real(sum(w_dsa.*cb,2))); I_dsa_f3(:,jf)=single(real(sum((w_dsa.*ct).*sgn_B300,2)));
    I_sa_t(:,jf)=single(real(sum(w_sa.*ct,2))); I_sa_f3(:,jf)=single(real(sum((w_sa.*ct).*sgn_B300,2)));
end
clear Cs_b1_all

iv_raw_b1 = core.band_inv_snr(I_raw_t, I_raw_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_raw_b3 = core.band_inv_snr(I_raw_t, I_raw_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b1 = core.band_inv_snr(I_dsa_t, I_dsa_f1, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_dsa_b3 = core.band_inv_snr(I_dsa_t, I_dsa_f3, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_sa_b3  = core.band_inv_snr(I_sa_t,  I_sa_f3,  freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_tr_raw = true_inv_snr(I_raw_s, I_raw_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);
iv_tr_dsa = true_inv_snr(I_dsa_s, I_dsa_n, freqs, P.nzi, P.nxi, P.dx, cfg.snr_band_hz, cfg.snr_box_m, cfg.snr_cap);

fa = (freqs >= 10.0 & freqs <= 35.0);
l2e = @(e, t) median(arrayfun(@(j) median(log2(max(double(e(w_deep_mask,j)),1e-6)./max(double(t(w_deep_mask,j)),1e-6))), find(fa)));
noise_acc = struct('RAW_B1_log2err', l2e(iv_raw_b1, iv_tr_raw), 'RAW_B300_log2err', l2e(iv_raw_b3, iv_tr_raw), ...
                   'DSA_B1_log2err', l2e(iv_dsa_b1, iv_tr_dsa), 'DSA_B300_log2err', l2e(iv_dsa_b3, iv_tr_dsa));
fprintf('=== Mid-SNR %sb (snr_db=%.1f dB) | Noise Acc: RAW B1=%+.2f B300=%+.2f | DSA B1=%+.2f B300=%+.2f ===\n', ...
    mname, snr_db, noise_acc.RAW_B1_log2err, noise_acc.RAW_B300_log2err, noise_acc.DSA_B1_log2err, noise_acc.DSA_B300_log2err);

mnames = {'RAW_R0','DSA_R0','RAW_UQW1_25','RAW_UQW1b_25','RAW_UQW1b_100','SA_UQW1b_25','DSA_PH_25','DSA_IQ70_25','DSA_W1b_25','DSA_UQW1_25','DSA_UQW1b_25'};
IS = struct(); IN = struct();
for k = 1:numel(mnames), nc = 1+(nm-1)*~endsWith(mnames{k},'_R0'); IS.(mnames{k})=zeros(nin,nc); IN.(mnames{k})=zeros(nin,nc); end
xmax = 1.05 * max(max(TS_flat(:)), max(TS_dsa(:))) / min(ms); xg = linspace(0, xmax, 1201); dxg = xg(2);

for jf = 1:nf
    f = freqs(jf); Cs_s = double(Cs_sig_all(:,:,jf)); Cs_n = double(Cs_nu_all(:,:,jf));
    T0_25 = uq_tab(f, xg, 0, 25); T1_25 = uq_tab(f, xg, sig, 25); T1_100 = uq_tab(f, xg, sig, 100);
    IS.RAW_R0 = IS.RAW_R0 + real(sum(Cs_s,2)); IN.RAW_R0 = IN.RAW_R0 + real(sum(Cs_n,2));
    IS.DSA_R0 = IS.DSA_R0 + real(sum(w_dsa.*Cs_s,2)); IN.DSA_R0 = IN.DSA_R0 + real(sum(w_dsa.*Cs_n,2));
    ivr1=double(iv_raw_b1(:,jf)); ivr3=double(iv_raw_b3(:,jf)); ivd1=double(iv_dsa_b1(:,jf)); ivd3=double(iv_dsa_b3(:,jf)); ivs3=double(iv_sa_b3(:,jf));
    for im = 1:nm
        u=(TS_flat/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A1=it(T1_25.A); C1=it(T1_25.C); P1=it(T1_25.P);
        H1=(A1+s70*P1.*C1+P1.*ivr1)./(1+s70*C1+ivr1); H3=(A1+s70*P1.*C1+P1.*ivr3)./(1+s70*C1+ivr3); Hs=(A1+s70*P1.*C1+P1.*ivs3)./(1+s70*C1+ivs3);
        IS.RAW_UQW1_25(:,im)=IS.RAW_UQW1_25(:,im)+real(sum(H1.*Cs_s,2)); IN.RAW_UQW1_25(:,im)=IN.RAW_UQW1_25(:,im)+real(sum(H1.*Cs_n,2));
        IS.RAW_UQW1b_25(:,im)=IS.RAW_UQW1b_25(:,im)+real(sum(H3.*Cs_s,2)); IN.RAW_UQW1b_25(:,im)=IN.RAW_UQW1b_25(:,im)+real(sum(H3.*Cs_n,2));
        IS.SA_UQW1b_25(:,im)=IS.SA_UQW1b_25(:,im)+real(sum((w_sa.*Hs).*Cs_s,2)); IN.SA_UQW1b_25(:,im)=IN.SA_UQW1b_25(:,im)+real(sum((w_sa.*Hs).*Cs_n,2));
        A1h=it(T1_100.A); C1h=it(T1_100.C); P1h=it(T1_100.P); Hh=(A1h+s70*P1h.*C1h+P1h.*ivr3)./(1+s70*C1h+ivr3);
        IS.RAW_UQW1b_100(:,im)=IS.RAW_UQW1b_100(:,im)+real(sum(Hh.*Cs_s,2)); IN.RAW_UQW1b_100(:,im)=IN.RAW_UQW1b_100(:,im)+real(sum(Hh.*Cs_n,2));

        u=(TS_dsa/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A0=it(T0_25.A); C0=it(T0_25.C); P0=it(T0_25.P); A1=it(T1_25.A); C1=it(T1_25.C); P1=it(T1_25.P);
        Hiq=(A0+s70*P0.*C0)./(1+s70*C0); Hw1b=(A0+s70*P0.*C0+P0.*ivd3)./(1+s70*C0+ivd3);
        Huq1=(A1+s70*P1.*C1+P1.*ivd1)./(1+s70*C1+ivd1); Huqb=(A1+s70*P1.*C1+P1.*ivd3)./(1+s70*C1+ivd3);
        IS.DSA_PH_25(:,im)=IS.DSA_PH_25(:,im)+real(sum((w_dsa.*P0).*Cs_s,2)); IN.DSA_PH_25(:,im)=IN.DSA_PH_25(:,im)+real(sum((w_dsa.*P0).*Cs_n,2));
        IS.DSA_IQ70_25(:,im)=IS.DSA_IQ70_25(:,im)+real(sum((w_dsa.*Hiq).*Cs_s,2)); IN.DSA_IQ70_25(:,im)=IN.DSA_IQ70_25(:,im)+real(sum((w_dsa.*Hiq).*Cs_n,2));
        IS.DSA_W1b_25(:,im)=IS.DSA_W1b_25(:,im)+real(sum((w_dsa.*Hw1b).*Cs_s,2)); IN.DSA_W1b_25(:,im)=IN.DSA_W1b_25(:,im)+real(sum((w_dsa.*Hw1b).*Cs_n,2));
        IS.DSA_UQW1_25(:,im)=IS.DSA_UQW1_25(:,im)+real(sum((w_dsa.*Huq1).*Cs_s,2)); IN.DSA_UQW1_25(:,im)=IN.DSA_UQW1_25(:,im)+real(sum((w_dsa.*Huq1).*Cs_n,2));
        IS.DSA_UQW1b_25(:,im)=IS.DSA_UQW1b_25(:,im)+real(sum((w_dsa.*Huqb).*Cs_s,2)); IN.DSA_UQW1b_25(:,im)=IN.DSA_UQW1b_25(:,im)+real(sum((w_dsa.*Huqb).*Cs_n,2));
    end
end

ref_2d = core.prep_image(I_ref, P.nzi, P.nxi, z_in, P.dx, cfg.bg_sigma_m, cfg.z_mute_m);
ref_hf = core.highpass_kz(ref_2d, P.dx, 10.0);
W_deep = {[180 1420], [560 950]};
rms_w = @(A, W) sqrt(mean(A(z_in>=W{2}(1)&z_in<=W{2}(2), x_in>=W{1}(1)&x_in<=W{1}(2)).^2, 'all'));
ref_rcf = log10(rms_w(ref_2d, W_targ) / max(rms_w(ref_2d, W_flank), 1e-30));
iz_w = find(z_in>=560 & z_in<=950); ix_w = find(x_in>=180 & x_in<=1420);

res = struct('model', [mname 'b'], 'seed', seed, 'snr_db', snr_db, 'm', ms, 'dip_stats', dip_stats, 'noise_acc', noise_acc, 'W_deep', struct(), 'elapsed_s', 0);
for k = 1:numel(mnames)
    key = mnames{k}; Sm = IS.(key); Nm = IN.(key); nc = size(Sm,2);
    e1=zeros(1,nc); e1h=e1; eamp=e1; snr=e1;
    for c=1:nc
        S2=core.prep_image(Sm(:,c),P.nzi,P.nxi,z_in,P.dx,cfg.bg_sigma_m,cfg.z_mute_m);
        N2=core.prep_image(Nm(:,c),P.nzi,P.nxi,z_in,P.dx,cfg.bg_sigma_m,cfg.z_mute_m);
        T2=S2+N2; Th=core.highpass_kz(T2,P.dx,10.0);
        e1(c)=core.score_e1(T2,ref_2d,iz_w,ix_w); e1h(c)=core.score_e1(Th,ref_hf,iz_w,ix_w);
        snr(c)=10*log10(sum(S2(iz_w,ix_w).^2,'all')/max(sum(N2(iz_w,ix_w).^2,'all'),1e-30));
        eamp(c)=abs(log10(rms_w(T2,W_targ)/max(rms_w(T2,W_flank),1e-30))-ref_rcf);
    end
    res.W_deep.(key)=struct('E1',e1,'E1_mean',mean(e1),'E1_HF',e1h,'E1_HF_mean',mean(e1h),'E_amp',eamp,'E_amp_mean',mean(eamp),'SNR_img_dB',snr,'SNR_img_dB_mean',mean(snr));
    fprintf('  %-15s | E1=%.4f (%.3f/%.3f/%.3f) | HF=%.4f | Eamp=%.3f | SNR=%+.2f dB\n', key, mean(e1), e1(1), e1(min(3,nc)), e1(end), mean(e1h), mean(eamp), mean(snr));
end
res.elapsed_s = toc(t0_all);
end

function res = run_u4b(core, cfg)
t0_all = tic; seed = 20261304; rng(seed, 'twister');
P.dx=10; P.npml=18; P.nzi=100; P.nxi=160; P.nz=136; P.nx=196;
P.z=((1:P.nz)'-19)*10; P.x=((1:P.nx)-19)*10; P.in={19:118, 19:178}; P.izs=21;
P.ixs=19+(8:3:152); P.ixr=19+(0:1:159); ns=numel(P.ixs); nr=numel(P.ixr); N=P.nz*P.nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2}); lin=sub2ind([P.nz P.nx],Zi(:),Xi(:)); nin=numel(lin);
isrc=sub2ind([P.nz P.nx],P.izs*ones(1,ns),P.ixs); irec=sub2ind([P.nz P.nx],P.izs*ones(1,nr),P.ixr);
[Zfull,Xfull]=ndgrid(P.z,P.x); z_in=P.z(P.in{1}); x_in=P.x(P.in{2}); [Zin,Xin]=ndgrid(z_in,x_in);
xs=P.x(P.ixs); xr=P.x(P.ixr);

rec_mask = false(nr, ns); xr_lo = zeros(1,ns); xr_hi = zeros(1,ns);
for s=1:ns, rec_mask(:,s)=(xr>=xs(s)-240)&(xr<=xs(s)+860); xr_lo(s)=min(xr(rec_mask(:,s))); xr_hi(s)=max(xr(rec_mask(:,s))); end
v_bg=1800+1.45*max(Zfull,0); z_dip=Zfull+0.11*(Xfull-800); dv=zeros(size(Zfull));
ze=[175,225,285,345,415,485,555,590,630,665,705,745,785,825,860,900,945];
for il=1:numel(ze)-1, dv(z_dip>=ze(il)&z_dip<ze(il+1))=max(min(0.055*randn(),0.09),-0.09); end
v_true=v_bg.*(1+dv); v_mig=core.smooth2d(v_true,3.0);
Q_true=55*ones(P.nz,P.nx); Q_true(Zfull<80)=1000; Q_true(Zfull>=520)=35;

s_mig=1./v_mig; eps_q=0.5; F.T=zeros(P.nz,P.nx,ns,'single'); F.ts=F.T;
for k=1:ns
    F.T(:,:,k)=eikonal_fs(s_mig,10,P.izs,P.ixs(k));
    F.ts(:,:,k)=(eikonal_fs(s_mig.*(1+eps_q./Q_true),10,P.izs,P.ixs(k))-eikonal_fs(s_mig.*(1-eps_q./Q_true),10,P.izs,P.ixs(k)))/(2*eps_q);
end
TS_flat=zeros(nin,ns); Fts_in=zeros(nin,ns);
for k=1:ns, [~,tt]=two_way(F,k,P); TS_flat(:,k)=double(tt(lin)); tk=F.ts(:,:,k); Fts_in(:,k)=double(tk(lin)); end; clear F

freqs=4:1:37; nf=numel(freqs); fref=25; fpeak_true=20; sig=cfg.sig; s70=gain_s(70); ms=[0.5 0.7 1.0 1.4 2.0]; nm=numel(ms);
corr_noise_sig=5.0; snr_db=2.0; rc=ceil(3*corr_noise_sig); g_nc=exp(-0.5*((-rc:rc)'/corr_noise_sig).^2); g_nc=g_nc/sqrt(sum(g_nc.^2));
rng(seed+1,'twister'); EPS_B1=2*(rand(nr,ns)>.5)-1; rng(seed+300,'twister'); sgn_B300=2*(rand(1,ns)>.5)-1;
W_mig=ones(1,nf); a1=freqs<6; W_mig(a1)=0.5*(1-cos(pi*(freqs(a1)-4)/2)); a2=freqs>34; W_mig(a2)=0.5*(1+cos(pi*(freqs(a2)-34)/3));

W20=(20/fpeak_true)^2*exp(-(20/fpeak_true)^2); b20=zeros(N,ns); b20(sub2ind([N ns],isrc,1:ns))=-W20/100;
[Lm,Um,Pm,Qm,Rm]=lu(helm_op(v_mig,20,P)); Ub=Qm*(Um\(Lm\(Pm*(Rm\b20)))); Dbg20=Ub(irec,:); clear Ub Lm Um Pm Qm Rm
vc20=v_true./(1-(1./(pi*Q_true))*log(20/fref)-1i./(2*Q_true));
[Lv,Uv,Pv,Qv,Rv]=lu(helm_op(vc20,20,P)); Uv20=Qv*(Uv\(Lv\(Pv*(Rv\b20)))); Dcal=(Uv20(irec,:)-Dbg20).*rec_mask; clear Uv20 Lv Uv Pv Qv Rv
sig_field=sqrt(mean(abs(Dcal(rec_mask)).^2))*10^(-snr_db/20);

Cs_sig_all=zeros(nin,ns,nf,'single'); Cs_nu_all=Cs_sig_all; Cs_b1_all=Cs_sig_all; I_RAW_R0=zeros(nin,1);
for jf=1:nf
    f=freqs(jf); w=2*pi*f; W_true=(f/fpeak_true)^2*exp(-(f/fpeak_true)^2);
    bt=zeros(N,ns); bt(sub2ind([N ns],isrc,1:ns))=-W_true/100;
    bm=zeros(N,ns); bm(sub2ind([N ns],isrc,1:ns))=-W_mig(jf)/100;
    [Lm,Um,Pm,Qm,Rm]=lu(helm_op(v_mig,f,P));
    Ub=Qm*(Um\(Lm\(Pm*(Rm\bt)))); D_bg=Ub(irec,:); clear Ub
    Um_=Qm*(Um\(Lm\(Pm*(Rm\bm)))); Us_c=conj(Um_(lin,:))*w^2; clear Um_
    vc=v_true./(1-(1./(pi*Q_true))*log(max(f,1e-6)/fref)-1i./(2*Q_true));
    [Lv,Uv,Pv,Qv,Rv]=lu(helm_op(vc,f,P)); Uv_=Qv*(Uv\(Lv\(Pv*(Rv\bt)))); D_vi=(Uv_(irec,:)-D_bg).*rec_mask; clear Uv_ Lv Uv Pv Qv Rv
    N_f=conv2(g_nc,1,sig_field*(1/sqrt(2))*(randn(nr,ns)+1i*randn(nr,ns)),'same').*rec_mask;
    Vs=extract_bp(D_vi,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10); Vn=extract_bp(N_f,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10);
    Vb=extract_bp((D_vi+N_f).*EPS_B1,Lm,Um,Pm,Qm,Rm,irec,lin,N,ns,10); clear Lm Um Pm Qm Rm
    cs=Us_c.*Vs; cn=Us_c.*Vn; cb=Us_c.*Vb;
    Cs_sig_all(:,:,jf)=single(cs); Cs_nu_all(:,:,jf)=single(cn); Cs_b1_all(:,:,jf)=single(cb);
    I_RAW_R0=I_RAW_R0+real(sum(cs+cn,2));
end

R0_2d=core.prep_image(I_RAW_R0,P.nzi,P.nxi,z_in,10,cfg.bg_sigma_m,cfg.z_mute_m);
[phi_eff,~]=core.blind_dip(R0_2d,10,cfg.dip_sigma_m,cfg.dip_cap_deg);
[w_dsa,TS_dsa]=core.build_dsa_in(Xin(:),Zin(:),20,xs,xr_lo,xr_hi,phi_eff(:),Fts_in,cfg);

I_raw_s=zeros(nin,nf,'single'); I_raw_n=I_raw_s; I_raw_t=I_raw_s; I_raw_f1=I_raw_s; I_raw_f3=I_raw_s;
I_dsa_s=I_raw_s; I_dsa_n=I_raw_s; I_dsa_t=I_raw_s; I_dsa_f1=I_raw_s; I_dsa_f3=I_raw_s;
for jf=1:nf
    cs=double(Cs_sig_all(:,:,jf)); cn=double(Cs_nu_all(:,:,jf)); ct=cs+cn; cb=double(Cs_b1_all(:,:,jf));
    I_raw_s(:,jf)=single(real(sum(cs,2))); I_raw_n(:,jf)=single(real(sum(cn,2))); I_raw_t(:,jf)=single(real(sum(ct,2)));
    I_raw_f1(:,jf)=single(real(sum(cb,2))); I_raw_f3(:,jf)=single(real(sum(ct.*sgn_B300,2)));
    I_dsa_s(:,jf)=single(real(sum(w_dsa.*cs,2))); I_dsa_n(:,jf)=single(real(sum(w_dsa.*cn,2))); I_dsa_t(:,jf)=single(real(sum(w_dsa.*ct,2)));
    I_dsa_f1(:,jf)=single(real(sum(w_dsa.*cb,2))); I_dsa_f3(:,jf)=single(real(sum((w_dsa.*ct).*sgn_B300,2)));
end; clear Cs_b1_all
iv_raw_b1=core.band_inv_snr(I_raw_t,I_raw_f1,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_raw_b3=core.band_inv_snr(I_raw_t,I_raw_f3,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_dsa_b1=core.band_inv_snr(I_dsa_t,I_dsa_f1,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_dsa_b3=core.band_inv_snr(I_dsa_t,I_dsa_f3,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_tr_raw=true_inv_snr(I_raw_s,I_raw_n,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
iv_tr_dsa=true_inv_snr(I_dsa_s,I_dsa_n,freqs,P.nzi,P.nxi,10,cfg.snr_band_hz,cfg.snr_box_m,cfg.snr_cap);
w_res_mask=(Zin(:)>=540&Zin(:)<=940)&(Xin(:)>=240&Xin(:)<=1360);
fa=(freqs>=10&freqs<=32);
l2e=@(e,t) median(arrayfun(@(j) median(log2(max(double(e(w_res_mask,j)),1e-6)./max(double(t(w_res_mask,j)),1e-6))), find(fa)));
noise_acc=struct('RAW_B1_log2err',l2e(iv_raw_b1,iv_tr_raw),'RAW_B300_log2err',l2e(iv_raw_b3,iv_tr_raw), ...
                 'DSA_B1_log2err',l2e(iv_dsa_b1,iv_tr_dsa),'DSA_B300_log2err',l2e(iv_dsa_b3,iv_tr_dsa));

anc_mask=(Zin(:)>=190&Zin(:)<=320)&(Xin(:)>=240&Xin(:)<=1360);
S_anc_raw=arrayfun(@(j) rms(double(I_raw_t(anc_mask,j))),1:nf);
S_anc_dsa=arrayfun(@(j) rms(double(I_dsa_t(anc_mask,j))),1:nf);
[W_eq_raw,beta_r,~]=core.fit_psae(freqs,S_anc_raw,[6 32],20);
[W_eq_dsa,beta_d,~]=core.fit_psae(freqs,S_anc_dsa,[6 32],20);
W_env_dsa=core.fit_env_anc(freqs,S_anc_dsa,W_mig,6.0,20);
tanc_d=median(TS_dsa(anc_mask,:),'all'); tanc_f=median(TS_flat(anc_mask,:),'all');
TS_dr=max(TS_dsa-tanc_d,0); TS_fr=max(TS_flat-tanc_f,0);

Img=struct('RAW_R0',zeros(nin,1),'DSA_R0',zeros(nin,1),'DSA_R0_Eq',zeros(nin,1), ...
           'RAW_UQW1_25',zeros(nin,nm),'RAW_UQW1b_25',zeros(nin,nm),'DSA_UQW1b_Unanc_25',zeros(nin,nm), ...
           'RAW_UQW1b_PAnc_25',zeros(nin,nm),'DSA_IQ70_PAnc_25',zeros(nin,nm),'DSA_W1b_PAnc_25',zeros(nin,nm), ...
           'DSA_UQW1_PAnc_25',zeros(nin,nm),'DSA_UQW1b_PAnc_25',zeros(nin,nm),'DSA_UQW1b_EnvAnc_25',zeros(nin,nm));
xmax=1.05*max(max(TS_flat(:)),max(TS_dsa(:)))/min(ms); xg=linspace(0,xmax,1201); dxg=xg(2);
for jf=1:nf
    f=freqs(jf); Cs_t=double(Cs_sig_all(:,:,jf))+double(Cs_nu_all(:,:,jf));
    wr=W_eq_raw(jf); wd=W_eq_dsa(jf); we=W_env_dsa(jf);
    raw_f=real(sum(Cs_t,2)); dsa_f=real(sum(w_dsa.*Cs_t,2));
    Img.RAW_R0=Img.RAW_R0+raw_f; Img.DSA_R0=Img.DSA_R0+dsa_f; Img.DSA_R0_Eq=Img.DSA_R0_Eq+wd*dsa_f;
    T0=uq_tab(f,xg,0,fref); T1=uq_tab(f,xg,sig,fref);
    ivr1=double(iv_raw_b1(:,jf)); ivr3=double(iv_raw_b3(:,jf)); ivd1=double(iv_dsa_b1(:,jf)); ivd3=double(iv_dsa_b3(:,jf));
    for im=1:nm
        u=(TS_flat/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A1=it(T1.A); C1=it(T1.C); P1=it(T1.P);
        Img.RAW_UQW1_25(:,im)=Img.RAW_UQW1_25(:,im)+real(sum(((A1+s70*P1.*C1+P1.*ivr1)./(1+s70*C1+ivr1)).*Cs_t,2));
        Img.RAW_UQW1b_25(:,im)=Img.RAW_UQW1b_25(:,im)+real(sum(((A1+s70*P1.*C1+P1.*ivr3)./(1+s70*C1+ivr3)).*Cs_t,2));
        u=(TS_dsa/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A1=it(T1.A); C1=it(T1.C); P1=it(T1.P);
        Img.DSA_UQW1b_Unanc_25(:,im)=Img.DSA_UQW1b_Unanc_25(:,im)+real(sum((w_dsa.*((A1+s70*P1.*C1+P1.*ivd3)./(1+s70*C1+ivd3))).*Cs_t,2));
        u=(TS_fr/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A1=it(T1.A); C1=it(T1.C); P1=it(T1.P);
        Img.RAW_UQW1b_PAnc_25(:,im)=Img.RAW_UQW1b_PAnc_25(:,im)+wr*real(sum(((A1+s70*P1.*C1+P1.*ivr3)./(1+s70*C1+ivr3)).*Cs_t,2));
        u=(TS_dr/ms(im))/dxg; i0=min(floor(u),numel(xg)-2); ww=u-i0; i0=i0+1; it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
        A0=it(T0.A); C0=it(T0.C); P0=it(T0.P); A1=it(T1.A); C1=it(T1.C); P1=it(T1.P);
        fiq=real(sum((w_dsa.*((A0+s70*P0.*C0)./(1+s70*C0))).*Cs_t,2));
        fw1=real(sum((w_dsa.*((A0+s70*P0.*C0+P0.*ivd3)./(1+s70*C0+ivd3))).*Cs_t,2));
        fu1=real(sum((w_dsa.*((A1+s70*P1.*C1+P1.*ivd1)./(1+s70*C1+ivd1))).*Cs_t,2));
        fub=real(sum((w_dsa.*((A1+s70*P1.*C1+P1.*ivd3)./(1+s70*C1+ivd3))).*Cs_t,2));
        Img.DSA_IQ70_PAnc_25(:,im)=Img.DSA_IQ70_PAnc_25(:,im)+wd*fiq;
        Img.DSA_W1b_PAnc_25(:,im)=Img.DSA_W1b_PAnc_25(:,im)+wd*fw1;
        Img.DSA_UQW1_PAnc_25(:,im)=Img.DSA_UQW1_PAnc_25(:,im)+wd*fu1;
        Img.DSA_UQW1b_PAnc_25(:,im)=Img.DSA_UQW1b_PAnc_25(:,im)+wd*fub;
        Img.DSA_UQW1b_EnvAnc_25(:,im)=Img.DSA_UQW1b_EnvAnc_25(:,im)+we*fub;
    end
end

well_xs=[400 600 800 1000 1200]; iz_res=find(z_in>=540&z_in<=940); v_in=v_true(P.in{1},P.in{2});
rw_map=zeros(size(v_in)); rw_map(2:end,:)=diff(log(v_in),1,1);
eval_w=@(v) score_w5(core, core.prep_image(v,P.nzi,P.nxi,z_in,10,cfg.bg_sigma_m,cfg.z_mute_m), rw_map, x_in, iz_res, well_xs);
res=struct('model','U4b','seed',seed,'snr_db',snr_db,'m',ms,'beta_anc_raw',beta_r,'beta_anc_dsa',beta_d,'noise_acc',noise_acc,'methods',struct(),'elapsed_s',0);
fprintf('=== Mid-SNR U4b Volve Twin | Noise Acc: RAW B1=%+.2f B300=%+.2f | DSA B1=%+.2f B300=%+.2f ===\n', ...
    noise_acc.RAW_B1_log2err, noise_acc.RAW_B300_log2err, noise_acc.DSA_B1_log2err, noise_acc.DSA_B300_log2err);
mkeys=fieldnames(Img)';
for ik=1:numel(mkeys)
    k=mkeys{ik}; M=Img.(k); nc=size(M,2); ts_=arrayfun(@(c) eval_w(M(:,c)), 1:nc);
    res.methods.(k)=struct('tie',ts_,'mean',mean(ts_));
    if nc==1, fprintf('  %-20s : %.4f\n', k, ts_(1));
    else, fprintf('  %-20s : mean=%.4f (%.3f/%.3f/%.3f)\n', k, mean(ts_), ts_(1), ts_(3), ts_(5)); end
end
res.elapsed_s=toc(t0_all);
end

function c_mean = score_w5(core, A, rw_map, x_in, iz_res, well_xs)
cs=zeros(1,numel(well_xs));
for iw=1:numel(well_xs)
    [~,jx]=min(abs(x_in-well_xs(iw))); tr=mean(A(iz_res,max(jx-1,1):min(jx+1,numel(x_in))),2);
    tb=core.bp_z(tr,10,4,18); rb=core.bp_z(rw_map(iz_res,jx),10,4,18);
    best=-Inf; for L=-3:3, best=max(best,core.cc(circshift(tb,L),rb)); end; cs(iw)=best;
end
c_mean=mean(cs);
end

function iv = true_inv_snr(Is_mat, In_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
[nin, nf] = size(Is_mat); bw = max(1, round(box_m / dx)); iv = zeros(nin, nf, 'single');
for j = 1:nf
    bnd = abs(freqs - freqs(j)) <= band_hz + 1e-9;
    ps = movmean(movmean(reshape(sum(double(Is_mat(:,bnd)),2),nzi,nxi).^2,bw,1),bw,2);
    pn = movmean(movmean(reshape(sum(double(In_mat(:,bnd)),2),nzi,nxi).^2,bw,1),bw,2);
    iv(:, j) = single(min(pn(:) ./ max(ps(:), realmin('single')), cap));
end
end

function V = extract_bp(Dmat, L, U, Pp, Qq, R, irec, lin, N, ns, dx)
b = zeros(N, ns); b(irec, :) = Dmat;
Vf = -(R \ (Pp' * (L' \ (U' \ (Qq' * b))))) / dx^2; V = Vf(lin, :);
end
