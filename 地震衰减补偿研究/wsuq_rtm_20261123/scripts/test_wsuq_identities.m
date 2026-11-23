function test_wsuq_identities()
% Unit self-test of ws_uq_apply.m verifying all 3 exact mathematical boundary identities (G0-1).
rng(20261123, 'twister');
nin = 200; ns = 51; sig = 0.447; f = 30; fref = 25; s70 = 1e-7;
Cs = randn(nin, ns) + 1i * randn(nin, ns);
iv = 0.05 + 2.0 * rand(nin, 1);

% Identity 1: TS == 0 => WS_UQ1 == R0
TS0 = zeros(nin, ns);
[R0, UQ0, WS0, WS1] = ws_uq_apply(Cs, TS0, f, sig, fref, iv, s70, []);
e1 = max(abs(WS1 - R0));
fprintf('ID1 (TS=0 -> WS_UQ1==R0): max|diff| = %.3e\n', e1);

% Identity 2: Equal TS across shots => WS_UQ1 == UQW1_old (fixed fref AND VA-UQ)
t_eq = repmat(linspace(0.002, 0.035, nin)', 1, ns);
[~, UQ_eq, ~, WS1_eq] = ws_uq_apply(Cs, t_eq, f, sig, fref, iv, s70, []);
e2_fix = max(abs(WS1_eq - UQ_eq)) / max(abs(UQ_eq));
[~, UQ_va, ~, WS1_va] = ws_uq_apply(Cs, t_eq, f, sig, fref, iv, s70, [15 35]);
e2_va = max(abs(WS1_va - UQ_va)) / max(abs(UQ_va));
fprintf('ID2 (equal t* -> WS_UQ1==UQW1): fixed=%.3e, VA=%.3e\n', e2_fix, e2_va);

% Identity 3: Heterogeneous TS (chimney undershooting t* in [11, 31] ms) => SNR gain > 5 dB
ts_row = 0.011 + 0.020 * exp(-(((0:ns-1) - 25)/10).^2);
TS_het = repmat(ts_row, nin, 1);
G_true = exp((-pi*f + 1i*2*f*log(f/fref)) * TS_het);
R_true = randn(nin, 1);
N_noise = 0.05 * (randn(nin, ns) + 1i * randn(nin, ns));
Cs_sig = G_true .* R_true;
Cs_tot = Cs_sig + N_noise;
iv_true = repmat(mean(abs(sum(N_noise,2)).^2) / mean(abs(sum(Cs_sig,2)).^2), nin, 1);
[~, UQ_h, WS0_h, WS1_h] = ws_uq_apply(Cs_tot, TS_het, f, sig, fref, iv_true, s70, []);
c_old = corr(UQ_h, R_true);
c_ws0 = corr(WS0_h, R_true);
c_ws1 = corr(WS1_h, R_true);
fprintf('ID3 (het t* 11-31ms, 30Hz): corr UQW1_old=%.4f, WS_UQ0=%.4f, WS_UQ1=%.4f\n', c_old, c_ws0, c_ws1);

ok = (e1 < 1e-6) && (e2_fix < 1e-5) && (e2_va < 1e-5) && (c_ws1 > c_old + 0.02);
if ~ok, error('test_wsuq_identities FAILED'); end
fprintf('ALL 3 BOUNDARY IDENTITIES PASSED!\n');
end
