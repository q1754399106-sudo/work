function core = uq_rtm_xmir_qrtm_core()
% UQ_RTM_XMIR_QRTM_CORE
% Upgraded heterogeneous bent-ray specular geometry (xmir_bent) +
% True Viscoacoustic Wave-Equation Q-RTM Helmholtz operators +
% Unified UQ-RTM & integrity-compliant evaluation helpers.
%
% Key Physics Upgrades over uq_rtm_unified_core.m:
%   1) build_dsa_bent: replaces homogeneous straight-ray xmir = x + dz*tan(th_s + 2*phi)
%      with exact 2D eikonal incident angle th_s = atan2(dT/dx, dT/dz) +
%      exact synchronized row-by-row 2nd-order Runge-Kutta (RK2) upward kinematic
%      ray tracing in heterogeneous v_mig(z,x), detecting turned/trapped rays and
%      computing exact surface emergence coordinate xmir_bent, surface emergence
%      angles, and bent-ray dip-consistent two-way attenuation TS_bent.
%   2) helm_op_visco: frequency-domain viscoacoustic Helmholtz operator with
%      Kjartansson/Kolsky causal complex velocity strictly inside the physical
%      interior and strictly real lossless C-PML absorbing strips (unconditionally
%      absorbing PML for both attenuating 'att' and anti-attenuating 'comp' modes).

core.default_cfg       = @default_cfg;
core.blind_dip         = @blind_dip;
core.build_dsa_str     = @build_dsa_str;
core.build_dsa_bent    = @build_dsa_bent;
core.trace_upward_rk2  = @trace_upward_rk2;
core.make_visco_vel    = @make_visco_vel;
core.helm_op_visco     = @helm_op_visco;
core.solve_helm_shot   = @solve_helm_shot;
core.solve_helm_adj    = @solve_helm_adj;
core.make_block_signs  = @make_block_signs;
core.band_inv_snr      = @band_inv_snr;
core.true_inv_snr      = @true_inv_snr;
core.fit_psae          = @fit_psae;
core.tab_idx           = @tab_idx;
core.prep_image        = @prep_image;
core.steep_smile_frac  = @steep_smile_frac;
core.res_band          = @res_band;
core.bp_z              = @bp_z;
core.load_well         = @load_well;
core.well_tie_full     = @well_tie_full;
core.cos_taper         = @cos_taper;
core.smooth2d          = @smooth2d;
core.highpass_kz       = @highpass_kz;
core.score_e1          = @score_e1;
core.cc                = @cc;
end

function cfg = default_cfg(preset)
if nargin < 1, preset = 'synth_10m'; end
cfg = struct();
cfg.B             = 300;
cfg.sig           = 0.447;
cfg.gamma         = 1.0;
cfg.fref          = 25.0;
cfg.dip_cap_deg   = 35.0;
cfg.snr_band_hz   = 5.0;
cfg.snr_cap       = 1e3;
cfg.w_floor       = 1e-4;
cfg.flip_seed     = 20261126;
switch preset
    case 'field_obc'
        cfg.anchor_mode   = 'psae';
        cfg.bg_sigma_m    = 40.0;
        cfg.z_mute_m      = 150.0;
        cfg.dip_sigma_m   = 60.0;
        cfg.arr_taper_m   = [-150.0, 250.0];
        cfg.theta_max_deg = 50.0;      % legacy straight-ray gate (unchanged for fair baseline)
        cfg.theta_tap_deg = 10.0;
        cfg.surf_max_deg  = 40.0;      % surface emergence angle cap (matches 42 deg Nyquist at 25 Hz, dxs=50m)
        cfg.inc_max_deg   = 55.0;      % local reflector incidence angle cap (mutes >55 deg overburden smiles)
        cfg.snr_box_m     = 400.0;
        cfg.anc_z_m       = [1000.0, 1600.0];
        cfg.anc_x_m       = [-2500.0, 1500.0];
        cfg.res_z_m       = [2600.0, 3200.0];
        cfg.res_x_m       = [-2500.0, 1500.0];
        cfg.ovb_z_m       = [1600.0, 2400.0];
    case 'synth_10m'
        cfg.anchor_mode   = 'none';
        cfg.bg_sigma_m    = 50.0;
        cfg.z_mute_m      = 160.0;
        cfg.dip_sigma_m   = 40.0;
        cfg.arr_taper_m   = [-80.0, 120.0];
        cfg.theta_max_deg = 75.0;
        cfg.theta_tap_deg = 10.0;
        cfg.surf_max_deg  = 70.0;
        cfg.inc_max_deg   = 75.0;
        cfg.snr_box_m     = 160.0;
        cfg.anc_z_m       = [190.0, 320.0];
        cfg.anc_x_m       = [240.0, 1360.0];
        cfg.res_z_m       = [560.0, 950.0];
        cfg.res_x_m       = [180.0, 1420.0];
    otherwise
        error('Unknown preset: %s', preset);
end
end

function [phi_eff, c_st, phi_raw] = blind_dip(R0_2d, dx, sig_m, cap_deg)
if nargin < 4, cap_deg = 35.0; end
sig_px = max(sig_m / dx, 1.0);
[gx, gz] = gradient(double(R0_2d), dx);
Jxx = smooth2d(gx.^2, sig_px);
Jzz = smooth2d(gz.^2, sig_px);
Jxz = smooth2d(gx.*gz, sig_px);
phi_raw = max(min(-0.5 * atan2(2*Jxz, Jzz - Jxx), deg2rad(cap_deg)), -deg2rad(cap_deg));
c_st    = sqrt((Jzz - Jxx).^2 + 4*Jxz.^2) ./ max(Jzz + Jxx, 1e-30);
phi_eff = c_st .* phi_raw;
end

function [w_dsa, TS_dsa, xmir_all, th_s_all, th_r_all] = build_dsa_str(xv, zv, z_src_vec, z_rec_val, xs, xr_lo, xr_hi, phi_v, F_ts_full, zi_v, xi_v, nz, nx, cfg)
% Legacy straight-ray homogeneous DSA geometry (exact reproduction of uq_rtm_unified_core.m)
nin = numel(xv); ns = numel(xs);
if numel(z_src_vec) == 1, z_src_vec = repmat(z_src_vec, ns, 1); end
[xs_sort, ord] = sort(xs(:));
w_dsa    = zeros(nin, ns, 'single');
TS_dsa   = zeros(nin, ns, 'single');
xmir_all = zeros(nin, ns, 'single');
th_s_all = zeros(nin, ns, 'single');
th_r_all = zeros(nin, ns, 'single');
dz_r     = max(zv(:) - z_rec_val, 10.0);
xv = xv(:); zv = zv(:); phi_v = phi_v(:);
zi_v = zi_v(:); xi_v = xi_v(:);
for k = 1:ns
    dz_s = max(zv - z_src_vec(k), 10.0);
    th_s = atan2(xv - xs(k), dz_s);
    th_r = max(min(th_s + 2*phi_v, deg2rad(78)), -deg2rad(78));
    xmir = xv + dz_r .* tan(th_r);
    d_in = min(xmir - xr_lo(k), xr_hi(k) - xmir);
    t_arr = cos_taper(d_in, cfg.arr_taper_m(1), cfg.arr_taper_m(2));
    t_gs  = cos_taper(cfg.theta_max_deg - abs(rad2deg(th_s)), 0.0, cfg.theta_tap_deg);
    t_gr  = cos_taper(cfg.theta_max_deg - abs(rad2deg(th_r)), 0.0, cfg.theta_tap_deg);
    w_dsa(:, k)    = single(max(t_arr .* t_gs .* t_gr, cfg.w_floor));
    xmir_all(:, k) = single(xmir);
    th_s_all(:, k) = single(rad2deg(th_s));
    th_r_all(:, k) = single(rad2deg(th_r));

    s0 = sub2ind([nz nx ns], zi_v, xi_v, k * ones(nin, 1));
    if ns == 1
        TS_dsa(:, k) = 2.0 * single(F_ts_full(s0));
    else
        q  = interp1(xs_sort, 1:ns, min(max(xmir, xs_sort(1)), xs_sort(end)));
        k0 = min(floor(q), ns - 1); ww = single(q - k0);
        a  = sub2ind([nz nx ns], zi_v, xi_v, ord(k0));
        b  = sub2ind([nz nx ns], zi_v, xi_v, ord(k0 + 1));
        TS_dsa(:, k) = single(F_ts_full(s0)) + (1 - ww).*single(F_ts_full(a)) + ww.*single(F_ts_full(b));
    end
end
end

function [w_bent, TS_bent, diag_info] = build_dsa_bent(v_mig, P, z_src_vec, z_rec_val, xs, xr_lo, xr_hi, phi_2d, F_T_full, F_ts_full, cfg)
% Upgraded 2D Heterogeneous Bent-Ray Specular DSA Geometry (xmir_bent):
%   1) Local bent-ray incident angle from eikonal gradient: th_s = atan2(dT/dx, dT/dz)
%      (smoothly blended to chord angle within r < 4*dx of source point).
%   2) Exact Snell reflection across blind dip phi(x): th_r_takeoff = th_s + 2*phi,
%      with top-side incidence check |th_s + phi| < 85 deg and upward check |th_r| < 85 deg.
%   3) Exact synchronized row-by-row 2nd-order Runge-Kutta (RK2) upward kinematic ray
%      tracing in heterogeneous slowness s(z,x) = 1/v_mig(z,x) from every interior node
%      to receiver depth z_rec_val, yielding exact xmir_bent, surface emergence angle,
%      and turning/trapped-ray flags.
%   4) Exact bent-ray dip-consistent two-way attenuation TS_bent = F_ts(src) + F_ts(xmir_bent).

nzi = P.nzi; nxi = P.nxi; nin = nzi * nxi; ns = numel(xs);
nz  = P.nz;  nx  = P.nx;  dx  = P.dx;
if numel(z_src_vec) == 1, z_src_vec = repmat(z_src_vec, ns, 1); end

z_in = P.z(P.in{1}); z_in = z_in(:);
x_in = P.x(P.in{2}); x_in = x_in(:)';
[Zin, Xin] = ndgrid(z_in, x_in);
[Zi_full, Xi_full] = ndgrid(P.in{1}, P.in{2});
zi_v = Zi_full(:); xi_v = Xi_full(:);

% Prepare full-grid slowness & lateral gradient for RK2 upward ray tracing
s_full = 1.0 ./ double(v_mig);
[sx_full, ~] = gradient(s_full, dx);
[~, iz_rec] = min(abs(P.z(:) - z_rec_val));
v_surf_row = double(v_mig(iz_rec, :));
x_full = P.x(:)';

[xs_sort, ord] = sort(xs(:));

w_bent     = zeros(nin, ns, 'single');
TS_bent    = zeros(nin, ns, 'single');
xmir_all   = zeros(nin, ns, 'single');
ths_all    = zeros(nin, ns, 'single');
thr_all    = zeros(nin, ns, 'single');
thrs_all   = zeros(nin, ns, 'single');
turned_all = false(nin, ns);

v_in = double(v_mig(P.in{1}, P.in{2}));

% For large field grids (nin*ns > 2e6), precompute shot-independent 3D RK2
% ray emergence table (z_c, x_c, theta_grid) ONCE (~5s) instead of re-tracing
% 124 shots x 220k pixels. Exact per-shot RK2 used on all smaller grids.
use_table = (double(nin) * double(ns) > 2.0e6);
if use_table
    sub = 2; % 20 m coarse grid (2.5x finer than 50 m vmig smoothing)
    iz_c = P.in{1}(1:sub:end); if iz_c(end) ~= P.in{1}(end), iz_c(end+1) = P.in{1}(end); end
    ix_c = P.in{2}(1:sub:end); if ix_c(end) ~= P.in{2}(end), ix_c(end+1) = P.in{2}(end); end
    th_grid = deg2rad(-82.5:2.5:82.5); n_th = numel(th_grid);
    nzc = numel(iz_c); nxc = numel(ix_c);
    vc_in = double(v_mig(iz_c, ix_c));
    X_tab  = zeros(nzc, nxc, n_th, 'single');
    PX_tab = zeros(nzc, nxc, n_th, 'single');
    TR_tab = zeros(nzc, nxc, n_th, 'single');
    for it = 1:n_th
        [xm_t, px_t, tr_t] = trace_upward_rk2( ...
            s_full, sx_full, iz_c, ix_c, iz_rec, x_full, dx, repmat(th_grid(it), nzc, nxc), vc_in);
        X_tab(:, :, it)  = single(xm_t);
        PX_tab(:, :, it) = single(px_t);
        TR_tab(:, :, it) = single(tr_t);
    end
    zc_vec = P.z(iz_c); xc_vec = P.x(ix_c);
    dth = th_grid(2) - th_grid(1); th0 = th_grid(1);
    % Bilinear spatial index maps from full interior (nzi x nxi) to coarse (nzc x nxc)
    rz_q = interp1(zc_vec, 1:nzc, z_in, 'linear', 'extrap');
    rx_q = interp1(xc_vec, 1:nxc, x_in, 'linear', 'extrap');
    rz_q = min(max(rz_q, 1), nzc - 1e-6); rx_q = min(max(rx_q, 1), nxc - 1e-6);
    iz0 = floor(rz_q(:)); wz = rz_q(:) - iz0; iz1 = min(iz0 + 1, nzc);
    ix0 = floor(rx_q(:)'); wx = rx_q(:)' - ix0; ix1 = min(ix0 + 1, nxc);
    % Pre-interpolate spatial slices to full interior? No, trilinear on the fly per shot is fast!
    [IZ0_m, IX0_m] = ndgrid(iz0, ix0); [IZ1_m, IX1_m] = ndgrid(iz1, ix1);
    [WZ_m,  WX_m ] = ndgrid(wz,  wx );
    w00 = single((1 - WZ_m) .* (1 - WX_m)); w10 = single(WZ_m .* (1 - WX_m));
    w01 = single((1 - WZ_m) .* WX_m);       w11 = single(WZ_m .* WX_m);
    b00 = sub2ind([nzc nxc], IZ0_m(:), IX0_m(:));
    b10 = sub2ind([nzc nxc], IZ1_m(:), IX0_m(:));
    b01 = sub2ind([nzc nxc], IZ0_m(:), IX1_m(:));
    b11 = sub2ind([nzc nxc], IZ1_m(:), IX1_m(:));
    stride_th = nzc * nxc;
end

for k = 1:ns
    % Step 1: Eikonal traveltime gradient -> exact bent-ray local incident angle
    Tk = double(F_T_full(:, :, k));
    [dTdx_f, dTdz_f] = gradient(Tk, dx);
    gx_T = dTdx_f(P.in{1}, P.in{2});
    gz_T = dTdz_f(P.in{1}, P.in{2});
    th_eik = atan2(gx_T, gz_T);

    % Blend to chord angle only within r < 4*dx of the source singularity
    dz_s = Zin - z_src_vec(k);
    dx_s = Xin - xs(k);
    r_s  = hypot(dx_s, dz_s);
    th_chord = atan2(dx_s, max(dz_s, 1.0));
    bl = cos_taper(r_s, 2.0*dx, 4.5*dx);
    th_s = bl .* th_eik + (1.0 - bl) .* th_chord;

    % Local reflector incidence & specular takeoff angle
    th_inc = th_s + phi_2d;
    th_r   = th_s + 2.0 * phi_2d;

    % Step 2: Synchronized 2D RK2 upward ray tracing (direct or via 3D emergence table)
    if ~use_table
        [xmir_2d, px_surf_2d, turned_2d] = trace_upward_rk2( ...
            s_full, sx_full, P.in{1}, P.in{2}, iz_rec, x_full, dx, th_r, v_in);
    else
        tr_bad = abs(th_r(:)) >= deg2rad(82.5);
        th_cl  = min(max(th_r(:), th0), th_grid(end) - 1e-6);
        u_th   = (th_cl - th0) / dth + 1.0;
        it0    = min(max(floor(u_th), 1), n_th - 1);
        wt     = single(u_th - it0); owt = single(1.0) - wt;
        off0   = (it0 - 1) * stride_th;
        off1   = it0 * stride_th;
        interp3_fast = @(T) owt .* (w00(:).*T(b00+off0) + w10(:).*T(b10+off0) + w01(:).*T(b01+off0) + w11(:).*T(b11+off0)) + ...
                            wt  .* (w00(:).*T(b00+off1) + w10(:).*T(b10+off1) + w01(:).*T(b01+off1) + w11(:).*T(b11+off1));
        xmir_2d    = reshape(double(interp3_fast(X_tab)),  nzi, nxi);
        px_surf_2d = reshape(double(interp3_fast(PX_tab)), nzi, nxi);
        turned_2d  = reshape(tr_bad | (interp3_fast(TR_tab) > 0.25), nzi, nxi);
    end

    % Surface emergence angles (governed by conserved horizontal slowness p_x)
    xmir_cl = min(max(xmir_2d, x_full(1)), x_full(end));
    v_at_mir = interp1(x_full, v_surf_row, xmir_cl, 'linear', 'extrap');
    th_r_surf = asin(min(max(px_surf_2d .* v_at_mir, -1.0), 1.0));
    % Source-side surface takeoff angle estimate via local p_x * v_src (exact in 1D v(z), accurate in 2D)
    [~, ix_sk] = min(abs(x_full - xs(k)));
    v_at_src = double(v_mig(min(max(round((z_src_vec(k) - P.z(1))/dx) + 1, 1), nz), ix_sk));
    px_s_loc = sin(th_s) ./ v_in;
    th_s_surf = asin(min(max(px_s_loc * v_at_src, -1.0), 1.0));

    % Step 3: Active-array + physical emergence + incidence/aperture soft gating
    d_in = min(xmir_2d - xr_lo(k), xr_hi(k) - xmir_2d);
    d_in(turned_2d) = -1e6;

    t_arr   = cos_taper(d_in, cfg.arr_taper_m(1), cfg.arr_taper_m(2));
    t_inc   = cos_taper(cfg.inc_max_deg - abs(rad2deg(th_inc)), 0.0, cfg.theta_tap_deg);
    t_up    = cos_taper(84.0 - abs(rad2deg(th_r)), 0.0, 8.0);
    t_s_srf = cos_taper(cfg.surf_max_deg - abs(rad2deg(th_s_surf)), 0.0, cfg.theta_tap_deg);
    t_r_srf = cos_taper(cfg.surf_max_deg - abs(rad2deg(th_r_surf)), 0.0, cfg.theta_tap_deg);

    wk = max(t_arr .* t_inc .* t_up .* t_s_srf .* t_r_srf, cfg.w_floor);
    w_bent(:, k)    = single(wk(:));
    xmir_all(:, k)  = single(xmir_2d(:));
    ths_all(:, k)   = single(rad2deg(th_s(:)));
    thr_all(:, k)   = single(rad2deg(th_r(:)));
    thrs_all(:, k)  = single(rad2deg(th_r_surf(:)));
    turned_all(:,k) = turned_2d(:);

    % Step 4: Exact bent-ray dip-consistent two-way t* from eikonal cube
    s0 = sub2ind([nz nx ns], zi_v, xi_v, k * ones(nin, 1));
    if ns == 1
        TS_bent(:, k) = 2.0 * single(F_ts_full(s0));
    else
        xm_v = xmir_2d(:);
        xm_v(turned_2d(:)) = xs(k); % safe clamp for zero-weight turned rays
        q  = interp1(xs_sort, 1:ns, min(max(xm_v, xs_sort(1)), xs_sort(end)));
        k0 = min(floor(q), ns - 1); ww = single(q - k0);
        a  = sub2ind([nz nx ns], zi_v, xi_v, ord(k0));
        b  = sub2ind([nz nx ns], zi_v, xi_v, ord(k0 + 1));
        TS_bent(:, k) = single(F_ts_full(s0)) + (1 - ww).*single(F_ts_full(a)) + ww.*single(F_ts_full(b));
    end
end

diag_info = struct('xmir', xmir_all, 'th_s_deg', ths_all, 'th_r_deg', thr_all, ...
                   'th_r_surf_deg', thrs_all, 'turned', turned_all);
end

function [xmir_2d, px_surf_2d, turned_2d] = trace_upward_rk2(s_full, sx_full, iz_in, ix_in, iz_rec, x_full, dx, th_r_2d, v_in)
% Synchronized row-by-row 2nd-order Runge-Kutta (RK2 / Midpoint) upward ray tracer
% in 2D heterogeneous slowness s(z, x) = 1 / v(z, x).
% Hamiltonian ODE integrated upward (z -> z - dz, dz = dx > 0):
%   dx/(-dz)  = +px / sqrt(s^2 - px^2)
%   dpx/(-dz) = +(s * ds/dx) / sqrt(s^2 - px^2)
% Exact for all 1D v(z) (conserves px exactly!) and 2nd-order accurate for 2D v(x,z).

nzi = numel(iz_in); nxi = numel(ix_in);
x0_full = x_full(1); nx_full = numel(x_full);
inv_dx  = 1.0 / dx;

% Precompute initial x, px, and non-upward flag on interior grid
x_init  = repmat(x_full(ix_in), nzi, 1);
px_init = sin(th_r_2d) ./ v_in;
bad0    = abs(th_r_2d) >= deg2rad(85.0);

% State arrays hold current (x, px) for all rows activated so far (rows iz_curr..nzi)
x_cur   = x_init;
px_cur  = px_init;
trn_cur = bad0;

% Minimum allowed cos(theta) before declaring ray turned/trapped
cos_min = cosd(86.5);

iz_bot = iz_in(end);
iz_top = iz_in(1);

% March upward row by row from bottom interior row to iz_rec
for iz_f = iz_bot:-1:(iz_rec + 1)
    % Interior rows active at depth index iz_f are those with iz_in >= iz_f (works for ANY stride!)
    r_start = find(iz_in >= iz_f, 1, 'first');
    if isempty(r_start), continue; end
    rows = r_start:nzi;

    xc  = x_cur(rows, :);
    pxc = px_cur(rows, :);
    tr  = trn_cur(rows, :);

    % Row slowness profiles at current z and midpoint z - dx/2
    s_r0  = s_full(iz_f, :);
    sx_r0 = sx_full(iz_f, :);
    s_rm  = 0.5 * (s_full(iz_f, :) + s_full(iz_f - 1, :));
    sx_rm = 0.5 * (sx_full(iz_f, :) + sx_full(iz_f - 1, :));

    % Stage 1 (k1 at current z)
    u1  = min(max((xc - x0_full) * inv_dx + 1.0, 1.0), nx_full - 1.0001);
    i1  = floor(u1); w1 = u1 - i1;
    s1  = s_r0(i1)  .* (1.0 - w1) + s_r0(i1 + 1)  .* w1;
    sx1 = sx_r0(i1) .* (1.0 - w1) + sx_r0(i1 + 1) .* w1;

    pz2_1 = s1.^2 - pxc.^2;
    tr = tr | (pz2_1 <= (s1 * cos_min).^2);
    pz1 = sqrt(max(pz2_1, (s1 * cos_min).^2));

    x_mid  = xc  + (0.5 * dx) * (pxc ./ pz1);
    px_mid = pxc + (0.5 * dx) * ((s1 .* sx1) ./ pz1);

    % Stage 2 (k2 at midpoint z - dx/2)
    u2  = min(max((x_mid - x0_full) * inv_dx + 1.0, 1.0), nx_full - 1.0001);
    i2  = floor(u2); w2 = u2 - i2;
    s2  = s_rm(i2)  .* (1.0 - w2) + s_rm(i2 + 1)  .* w2;
    sx2 = sx_rm(i2) .* (1.0 - w2) + sx_rm(i2 + 1) .* w2;

    pz2_2 = s2.^2 - px_mid.^2;
    tr = tr | (pz2_2 <= (s2 * cos_min).^2);
    pz2 = sqrt(max(pz2_2, (s2 * cos_min).^2));

    xn  = xc  + dx * (px_mid ./ pz2);
    pxn = pxc + dx * ((s2 .* sx2) ./ pz2);

    % Keep turned rays frozen cleanly
    xn(tr)  = xc(tr);
    pxn(tr) = pxc(tr);

    x_cur(rows, :)   = xn;
    px_cur(rows, :)  = pxn;
    trn_cur(rows, :) = tr;
end

% Rows above iz_rec (if any, e.g. water layer above seabed receivers): straight fallback
above = (iz_in <= iz_rec);
if any(above)
    x_cur(above, :)   = x_init(above, :);
    px_cur(above, :)  = px_init(above, :);
    trn_cur(above, :) = false;
end

xmir_2d    = x_cur;
px_surf_2d = px_cur;
turned_2d  = trn_cur;
end

function vc = make_visco_vel(v_real, Q_est, f, fref, npml, mode, w_amp)
% Build Kjartansson/Kolsky complex velocity STRICTLY inside the physical interior;
% PML border strips remain strictly lossless real v_real for unconditional PML stability.
if nargin < 7 || isempty(w_amp), w_amp = 1.0; end
[nz, nx] = size(v_real);
in_mask = false(nz, nx);
in_mask((npml+1):(nz-npml), (npml+1):(nx-npml)) = true;

invQ = zeros(nz, nx);
finQ = isfinite(Q_est) & (Q_est > 0);
invQ(finQ) = 1.0 ./ double(Q_est(finQ));

% Causal dispersion slowness shift (active everywhere finite Q, clamped for safety)
ds_disp = -(invQ / pi) * log(max(f, 1e-6) / fref);
re_fac  = min(max(1.0 + ds_disp, 0.5), 2.0);

% Attenuation imaginary slowness shift (strictly zero inside PML!)
ds_att = (0.5 * w_amp) * invQ;
ds_att(~in_mask) = 0.0;

switch mode
    case 'acoustic'
        vc = double(v_real);
    case 'disp'
        vc = double(v_real) ./ re_fac;
    case 'att'
        vc = double(v_real) ./ (re_fac - 1i * ds_att);
    case 'comp'
        vc = double(v_real) ./ (re_fac + 1i * ds_att);
    otherwise
        error('Unknown viscoacoustic mode: %s', mode);
end
end

function A = helm_op_visco(vc, v_real, f, P)
% Viscoacoustic Helmholtz operator A = Lx + Lz + (w / vc)^2 with strictly REAL
% absorbing C-PML damping coefficient s0 computed from real velocity v_real.
[nz, nx] = size(vc); h = P.dx; w = 2*pi*f; np = P.npml; L = np*h;
s0 = 3.0 * max(double(v_real(:))) * log(1e4) / (2.0 * L);
Az = lap1_visco(nz, np, h, w, s0, L);
Bx = lap1_visco(nx, np, h, w, s0, L);
A  = kron(speye(nx), Az) + kron(Bx, speye(nz)) + spdiags((w ./ vc(:)).^2, 0, nz*nx, nz*nx);
end

function D = lap1_visco(n, np, h, w, s0, L)
d = @(y) min(max(max(np+1-y, y-(n-np)), 0), np) * h;
s = @(y) 1.0 - 1i * s0 * (d(y)/L).^2 / w;
i = (1:n)'; si = s(i); sp = s(i+0.5); sm = s(i-0.5);
pml = d(i)>0 | d(i-1)>0 | d(i+1)>0 | d(i-2)>0 | d(i+2)>0;
c_m = 1.0 ./ (si .* sm) / h^2;
c_p = 1.0 ./ (si .* sp) / h^2;
c_0 = -(c_m + c_p);
e = ones(n, 1);
q = [-e/12, 4*e/3, -5*e/2, 4*e/3, -e/12] / h^2;
B = zeros(n, 5);
B(pml, 2) = c_m(pml); B(pml, 3) = c_0(pml); B(pml, 4) = c_p(pml);
B(~pml, :) = q(~pml, :);
r = []; c = []; vals = []; off = -2:2;
for k = 1:5
    j = i + off(k); ok = j >= 1 & j <= n & B(:, k) ~= 0;
    r = [r; i(ok)]; c = [c; j(ok)]; vals = [vals; B(ok, k)]; %#ok<AGROW>
end
D = sparse(r, c, vals, n, n);
end

function U = solve_helm_shot(L, Uf, Pp, Qq, R, b)
U = Qq * (Uf \ (L \ (Pp * (R \ b))));
end

function V = solve_helm_adj(L, Uf, Pp, Qq, R, b, dx)
V = -(R \ (Pp' * (L' \ (Uf' \ (Qq' * b))))) / dx^2;
end

function [EPS_tr, shot_sgn] = make_block_signs(shot_of, rec_ix, ns, B, flip_seed)
ntrace = numel(shot_of);
if B == 1
    rng(flip_seed, 'twister');
    EPS_tr = single(2*(rand(ntrace, 1) > 0.5) - 1);
    shot_sgn = [];
elseif B >= 300
    rng(flip_seed + B, 'twister');
    shot_sgn = single(2*(rand(ns, 1) > 0.5) - 1);
    EPS_tr = shot_sgn(shot_of(:));
else
    nper = accumarray(shot_of(:), 1, [ns 1]);
    nb = ceil(max(nper) / B);
    rng(flip_seed + B, 'twister');
    s = single(2*(rand(nb, ns) > 0.5) - 1);
    blk = zeros(ntrace, 1);
    for k = 1:ns
        t = find(shot_of == k);
        [~, o] = sort(rec_ix(t));
        r = zeros(numel(t), 1); r(o) = 1:numel(t);
        blk(t) = ceil(r / B);
    end
    EPS_tr = s(sub2ind([nb ns], blk, shot_of(:)));
    shot_sgn = [];
end
end

function iv = band_inv_snr(I_mat, If_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
if nargin < 9, cap = 1e3; end
[nin, nf] = size(I_mat);
bw = max(1, round(box_m / dx));
iv = zeros(nin, nf, 'single');
for j = 1:nf
    bnd = abs(freqs - freqs(j)) <= band_hz + 1e-9;
    s_tot = sum(double(I_mat(:, bnd)), 2);
    s_flp = sum(double(If_mat(:, bnd)), 2);
    p  = box_pow(s_tot, nzi, nxi, bw);
    pf = box_pow(s_flp, nzi, nxi, bw);
    iv(:, j) = single(min(pf ./ max(p - pf, realmin('single')), cap));
end
end

function iv = true_inv_snr(Is_mat, In_mat, freqs, nzi, nxi, dx, band_hz, box_m, cap)
if nargin < 9, cap = 1e3; end
[nin, nf] = size(Is_mat);
bw = max(1, round(box_m / dx));
iv = zeros(nin, nf, 'single');
for j = 1:nf
    bnd = abs(freqs - freqs(j)) <= band_hz + 1e-9;
    ps = box_pow(sum(double(Is_mat(:, bnd)), 2), nzi, nxi, bw);
    pn = box_pow(sum(double(In_mat(:, bnd)), 2), nzi, nxi, bw);
    iv(:, j) = single(min(pn ./ max(ps, realmin('single')), cap));
end
end

function y = box_pow(v, nzi, nxi, bw)
y = reshape(v, nzi, nxi).^2;
y = movmean(movmean(y, bw, 1), bw, 2);
y = y(:);
end

function [W_eq, beta_anc, eff_tstar_ms] = fit_psae(freqs, S_anc, f_fit_range, f_mid)
fm = (freqs >= f_fit_range(1)) & (freqs <= f_fit_range(2));
p = polyfit(freqs(fm), log(max(S_anc(fm), 1e-30)), 1);
beta_anc = p(1);
eff_tstar_ms = -beta_anc / pi * 1e3;
W_eq = exp(-beta_anc * (freqs - f_mid));
end

function [idx, ww] = tab_idx(ts_m, dxg, ntab_m2)
u = ts_m / dxg;
i0 = min(max(floor(u), 0), ntab_m2);
ww = u - i0;
idx = int32(i0) + 1;
end

function y = prep_image(v, nzi, nxi, z_in, dx, bg_sigma_m, z_mute_m)
y = reshape(double(v), nzi, nxi);
s = bg_sigma_m / dx;
m = ceil(3 * s);
h = exp(-0.5 * ((-m:m) / s).^2); h = h / sum(h);
yp = [repmat(y(1,:), m, 1); y; repmat(y(end,:), m, 1)];
y = y - conv2(h', 1, yp, 'valid');
y(z_in < z_mute_m, :) = 0;
end

function fr = steep_smile_frac(A_mid, dx, deg_cut, k_lo_cyckm, k_hi_cyckm)
if nargin < 3, deg_cut = 25.0; end
if nargin < 4, k_lo_cyckm = 3.0; end
if nargin < 5, k_hi_cyckm = 25.0; end
[nz, nx] = size(A_mid);
win = hann(nz) * hann(nx)';
F2 = abs(fftshift(fft2(A_mid .* win))).^2;
kz = ((-floor(nz/2)):(ceil(nz/2)-1))' / (nz * dx);
kx = ((-floor(nx/2)):(ceil(nx/2)-1))  / (nx * dx);
[KZ, KX] = ndgrid(kz, kx);
k_mag = sqrt(KX.^2 + KZ.^2);
in_band = (k_mag >= k_lo_cyckm / 1000) & (k_mag <= k_hi_cyckm / 1000);
steep = in_band & (abs(KX) > tand(deg_cut) * max(abs(KZ), 1e-9));
fr = sum(F2(steep), 'all') / max(sum(F2(in_band), 'all'), 1e-30);
end

function [bw10, cen] = res_band(A_res, dx, nfft)
if nargin < 3, nfft = 256; end
kz = (0:nfft/2) / (nfft * dx) * 1000;
s = mean(abs(fft(A_res .* hann(size(A_res, 1)), nfft, 1)), 2);
s = s(1:nfft/2+1)';
[mx, i] = max(s); t = mx * 10^(-10/20);
l = find(s(1:i) < t, 1, 'last'); if isempty(l), l = 1; end
h = i - 1 + find(s(i:end) < t, 1); if isempty(h), h = numel(s); end
bw10 = [kz(l), kz(h)];
cen = sum(kz .* (s.^2)) / max(sum(s.^2), 1e-30);
end

function y = bp_z(r, dz, lo, hi)
n = 2^nextpow2(4 * numel(r));
R = fft(r, n);
k = abs([0:n/2, -n/2+1:-1]') / (n * dz) * 1000;
y = real(ifft(R .* (k >= lo & k <= hi)));
y = y(1:numel(r));
end

function w = load_well(fn, z, x)
fid = fopen(fn); a = fread(fid, Inf, 'single'); fclose(fid);
a = reshape(a, 4, [])';
w.zw = double(a(:, 1)); w.sw = double(a(:, 2)); w.rw = double(a(:, 4));
w.dzw = median(diff(w.zw));
[~, w.iz] = min(abs(z(:) - w.zw'), [], 1);
w.ix = arrayfun(@(xv) find(abs(x - xv) == min(abs(x - xv)), 1), w.sw);
end

function [best_cc, best_lag_m, zero_cc, tr_b] = well_tie_full(A, w, rb, dx_grid, lo, hi, half_avg_m, max_lag_m)
if nargin < 7, half_avg_m = 20.0; end
if nargin < 8, max_lag_m  = 30.0; end
na = round(half_avg_m / dx_grid);
nw = numel(w.zw); tr = zeros(nw, 1);
for i = 1:nw
    tr(i) = mean(A(w.iz(i), max(w.ix(i) - na, 1):min(w.ix(i) + na, size(A, 2))));
end
tr_b = bp_z(tr, w.dzw, lo, hi);
zero_cc = cc(tr_b, rb);
ns_ = round(max_lag_m / w.dzw);
best_cc = -Inf; best_lag_m = 0;
for L = -ns_:ns_
    c = cc(circshift(tr_b, L), rb);
    if c > best_cc
        best_cc = c;
        best_lag_m = L * w.dzw;
    end
end
end

function t = cos_taper(x, x0, x1)
u = min(max((x - x0) / max(x1 - x0, 1e-12), 0), 1);
t = 0.5 - 0.5 * cos(pi * u);
end

function B = smooth2d(A, sig)
r = max(1, ceil(3 * sig));
g = exp(-0.5 * ((-r:r) / sig).^2); g = g / sum(g);
Ap = [repmat(A(1,:), r, 1); A; repmat(A(end,:), r, 1)];
B1 = conv2(g', 1, Ap, 'valid');
B1p = [repmat(B1(:,1), 1, r), B1, repmat(B1(:,end), 1, r)];
B = conv2(1, g, B1p, 'valid');
end

function yh = highpass_kz(y, dx, kz_cut)
n = 2^nextpow2(2 * size(y, 1));
Y = fft(y, n, 1);
kz = abs([0:n/2, -n/2+1:-1]') / (n * dx) * 1000;
g = 0.5 - 0.5 * cos(pi * min(max((kz - (kz_cut - 2)) / 4, 0), 1));
yh = real(ifft(Y .* g, n, 1));
yh = yh(1:size(y, 1), :);
end

function s = score_e1(a, r, iz, ix)
rr = r(iz, ix);
c0 = cc(a(iz, ix), rr);
cl = c0;
nz = size(a, 1);
for L = [-2 -1 1 2]
    cl = max(cl, cc(a(min(max(iz + L, 1), nz), ix), rr));
end
s = 0.5 * (c0 + cl);
end

function c = cc(a, b)
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
c = (a' * b) / max(sqrt((a' * a) * (b' * b)), 1e-30);
end
