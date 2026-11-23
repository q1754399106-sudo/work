function step0_tstar_audit()
% T2.0a: Zero-LU inter-shot t* dispersion and Theorem-4 Cauchy-Schwarz variance-reduction ratio audit
% on cached Volve F-5 (absolute + anchored 2.27 km) and synthetic D2, D4 (v_band) t* fields.
% Read-only access to main caches. Low RAM (< 500 MB), 2 threads.

maxNumCompThreads(2);
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
out_dir = fullfile(root, 'results');
if ~isfolder(out_dir), mkdir(out_dir); end

main_res = 'D:\ChatGPT\研究\地震衰减补偿研究';
addpath(fullfile(main_res, 'uq_rtm_20261117', 'scripts'));

res = struct();

% ---- 1. Volve F-5 (absolute t* and anchored 2.27 km Delta t*) ----
fld_dir = fullfile(main_res, 'uq_rtm_20261117', 'cache', 'fld');
meta = jsondecode(fileread(fullfile(fld_dir, 'meta.json')));
fid = fopen(fullfile(fld_dir, 'src_ix.bin'), 'r'); src_ix = fread(fid, inf, 'int32'); fclose(fid);
dx = double(meta.dx); npml = double(meta.npml);
nz = double(meta.nz); nx = double(meta.nx);
nzi = double(meta.nzi); nxi = double(meta.nxi);
z_all = ((1:nz)' - npml - 1) * dx;
x_all = double(meta.xmin) + ((1:nx) - npml - 1) * dx;
in_z = npml + (1:nzi);
in_x = npml + (1:nxi);
z_in = z_all(in_z); x_in = x_all(in_x);
src_x = x_all(double(src_ix));
[Zi, Xi] = ndgrid(z_in, x_in);

% Masks: F-5 well reservoir window (x in [-1300, -1100], z in [2600, 3200]) and full eval reservoir
m_f5  = (Zi >= 2600 & Zi <= 3200) & (Xi >= -1300 & Xi <= -1100);
m_res = (Zi >= 2600 & Zi <= 3200) & (Xi >= -2500 & Xi <= 1500);

S_abs = load(fullfile(fld_dir, 'ts.mat'), 'TS');
TS_abs = double(S_abs.TS); clear S_abs
S_anc = load(fullfile(main_res, 'va_uq_field2_20261121', 'cache', 'anc', 'ts.mat'), 'TS');
TS_anc = double(S_anc.TS); clear S_anc

% Verify shape matches numel(Zi) x numel(src_x)
fprintf('Volve TS_abs: %d x %d (grid %d x %d = %d, ns=%d)\n', size(TS_abs,1), size(TS_abs,2), numel(z_in), numel(x_in), numel(Zi), numel(src_x));
res.volve_abs_F5  = audit_ts(TS_abs, Zi(:), Xi(:), src_x(:)', m_f5(:),  [15 22 35 45], [1.0 0.5], 'Volve_abs_F5');
res.volve_abs_res = audit_ts(TS_abs, Zi(:), Xi(:), src_x(:)', m_res(:), [15 22 35 45], [1.0 0.5], 'Volve_abs_res');
res.volve_anc_F5  = audit_ts(TS_anc, Zi(:), Xi(:), src_x(:)', m_f5(:),  [15 22 35 45], [1.0 0.5], 'Volve_anc227_F5');
clear TS_abs TS_anc

% ---- 2. Synthetic D2 & D4 (va_uq_step1 v_band t* fields) ----
P = rtm_setup();
[Z2, X2] = ndgrid(P.z(P.in{1}), P.x(P.in{2}));
lin = sub2ind([P.nz P.nx], repmat(P.in{1}(:), 1, numel(P.in{2})), repmat(P.in{2}(:)', numel(P.in{1}), 1));
lin = lin(:);
xs_syn = P.x(P.ixs);
m_wdeep = (Z2 >= 820 & Z2 <= 1550) & (X2 >= 300 & X2 <= 2100);
m_chim  = (Z2 >= 820 & Z2 <= 1550) & (X2 >= 950 & X2 <= 1450);

for nm = {'D2', 'D4'}
    name = nm{1};
    ff = fullfile(main_res, 'va_uq_step1_20261121', 'cache', ['fields_' name '_vb.mat']);
    F = load(ff);
    ns = numel(P.ixs);
    TS = zeros(numel(lin), ns);
    for s = 1:ns
        [~, tt] = two_way(F, s, P);
        TS(:, s) = double(tt(lin));
    end
    clear F
    res.([name '_Wdeep']) = audit_ts(TS, Z2(:), X2(:), xs_syn(:)', m_wdeep(:), [20 35 45], [1.0 0.5], [name '_Wdeep']);
    if strcmp(name, 'D2')
        res.D2_chimney = audit_ts(TS, Z2(:), X2(:), xs_syn(:)', m_chim(:), [20 35 45], [1.0 0.5], 'D2_chimney');
    end
    clear TS
end

% Gate G0-2 verdict (PLAN_ZH.md section 2: active aperture |xs-x| <= 1.5*z)
rcs_d2_35_m1  = res.D2_chimney.active.f35_m10.Rcs_med;
rcs_vf5_22_m05 = res.volve_abs_F5.active.f22_m05.Rcs_med;
rcs_vf5_35_m1  = res.volve_abs_F5.active.f35_m10.Rcs_med;
g0_2_pass = (rcs_d2_35_m1 >= 1.30) || (rcs_vf5_22_m05 >= 1.30) || (rcs_vf5_35_m1 >= 1.30);
res.verdict_G0_2 = struct('pass', g0_2_pass, 'D2_chim_f35_m1', rcs_d2_35_m1, ...
                          'Volve_F5_f22_m05', rcs_vf5_22_m05, 'Volve_F5_f35_m1', rcs_vf5_35_m1, 'threshold', 1.30);

txt = jsonencode(res, 'PrettyPrint', true);
fid = fopen(fullfile(out_dir, 'step0_tstar_audit.json'), 'w'); fwrite(fid, txt, 'char'); fclose(fid);
fprintf('\n=== VERDICT G0-2 ===\n');
fprintf('D2_chimney (35 Hz, m=1): R_CS med = %.3f (%.2f dB)\n', rcs_d2_35_m1, 10*log10(rcs_d2_35_m1));
fprintf('Volve_F5   (22 Hz, m=0.5): R_CS med = %.3f (%.2f dB)\n', rcs_vf5_22_m05, 10*log10(rcs_vf5_22_m05));
fprintf('Volve_F5   (35 Hz, m=1): R_CS med = %.3f (%.2f dB)\n', rcs_vf5_35_m1, 10*log10(rcs_vf5_35_m1));
fprintf('G0-2 PASS = %d\n', g0_2_pass);
end

function out = audit_ts(TS, Z, X, src_x, mask, freqs, ms, label)
idx = find(mask);
% Subsample at most 2000 pixels evenly if huge (for speed; exact on F5)
if numel(idx) > 2000, idx = idx(round(linspace(1, numel(idx), 2000))); end
Tsub = TS(idx, :);
Zsub = Z(idx);
Xsub = X(idx);
np = numel(idx);

fprintf('\n--- %s (n_pix=%d, ns=%d) ---\n', label, np, size(TS,2));
for mode = {'all', 'active'}
    md = mode{1};
    stds = zeros(np, 1); ranges = zeros(np, 1); p10s = zeros(np, 1); p90s = zeros(np, 1); nact = zeros(np, 1);
    for i = 1:np
        if strcmp(md, 'all')
            row = Tsub(i, :);
        else
            k = abs(src_x - Xsub(i)) <= 1.5 * Zsub(i);
            row = Tsub(i, k);
        end
        nact(i) = numel(row);
        stds(i) = std(row);
        ranges(i) = max(row) - min(row);
        q = prctile(row, [10 90]);
        p10s(i) = q(1); p90s(i) = q(2);
    end
    st = struct('n_shots_med', median(nact), 'std_ms_med', 1e3*median(stds), 'range_ms_med', 1e3*median(ranges), ...
                'p10_ms_med', 1e3*median(p10s), 'p90_ms_med', 1e3*median(p90s));
    for f = freqs
        for m = ms
            rcs = zeros(np, 1);
            wratio_old = zeros(np, 1);
            wratio_ws  = zeros(np, 1);
            for i = 1:np
                if strcmp(md, 'all'), row = Tsub(i, :) / m;
                else, row = Tsub(i, abs(src_x - Xsub(i)) <= 1.5 * Zsub(i)) / m;
                end
                e_neg = exp(-2 * pi * f * row);
                e_pos = exp(+2 * pi * f * row);
                rcs(i) = mean(e_neg) * mean(e_pos);
                % Weight ratio of top-25% high-t* shots vs bottom-25% low-t* shots
                rs = sort(row);
                nq = max(1, floor(numel(rs)/4));
                lo = rs(1:nq); hi = rs(end-nq+1:end);
                wratio_old(i) = mean(exp(+pi * f * hi)) / mean(exp(+pi * f * lo));
                wratio_ws(i)  = mean(exp(-pi * f * hi)) / mean(exp(-pi * f * lo));
            end
            tag = sprintf('f%d_m%02d', round(f), round(10*m));
            st.(tag) = struct('Rcs_med', median(rcs), 'Rcs_p90', prctile(rcs, 90), ...
                              'dB_med', 10*log10(median(rcs)), 'dB_p90', 10*log10(prctile(rcs, 90)), ...
                              'hi_over_lo_old_med', median(wratio_old), 'hi_over_lo_ws_med', median(wratio_ws));
        end
    end
    out.(md) = st;
    fprintf('  [%6s] nshot=%.0f, std=%.2f ms, range=%.2f ms (P10=%.1f, P90=%.1f ms)\n', ...
            md, st.n_shots_med, st.std_ms_med, st.range_ms_med, st.p10_ms_med, st.p90_ms_med);
    for f = freqs
        t1 = sprintf('f%d_m10', round(f)); t05 = sprintf('f%d_m05', round(f));
        fprintf('     f=%2d Hz | m=1: R_CS=%.2f (%.1f dB, P90=%.1f dB), hi/lo old=%.2fx ws=%.2fx | m=0.5: R_CS=%.2f (%.1f dB)\n', ...
                f, st.(t1).Rcs_med, st.(t1).dB_med, st.(t1).dB_p90, st.(t1).hi_over_lo_old_med, st.(t1).hi_over_lo_ws_med, ...
                st.(t05).Rcs_med, st.(t05).dB_med);
    end
end
end
