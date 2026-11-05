%% ei_qcomp_net.m
% Self-supervised Q compensation with Equivariant Imaging (EI) -- network version
% Requires: MATLAB R2021b+ with Deep Learning Toolbox (convolution1dLayer, dlarray)
%
% Compares on synthetic data:
%   1) Tikhonov (mu tuned with ground truth -- oracle classical baseline)
%   2) MC            : data consistency only          ||A f(y) - y||^2
%   3) MC + EI       : + shift equivariance            ||f(A T_g x^) - T_g x^||^2
%   4) MC + EI + fmax: + physics-derived time-varying band limit f_max(t)
%
% Options below: Q error test (Qerr), noise level, iterations, EI weight.
% Environment variables (read at start): EI_SEED, EI_QERR, EI_ITERS, EI_CFG (e.g. "1 2"),
% EI_SCHED (cosine | const), EI_SMOKE=1 for a 30-step debug run,
% EI_SCEN (stationary | density | amp | broad2narrow | narrow2broad) for depth-varying reflectivity,
% EI_FREF (Hz) for the dispersion reference frequency of both data and operator (default Nyquist).
% Python/numpy linear prototype: ei_qcomp_linear.py

clear;

%% ---------------- parameters ----------------
P.N      = 256;      % samples per trace
P.dt     = 0.004;    % s
P.Q      = 50;       % true Q (used to generate data)
P.Qerr   = 0.0;      % relative Q error in the operator used for training, e.g. 0.2 -> Q_model = 1.2*Q
P.fp     = 25;       % Ricker peak frequency (Hz)
P.snrdB  = 30;       % noise level of attenuated data
P.ntr    = 4000;     % training traces (no labels used except for evaluation / baseline tuning)
P.nte    = 500;      % test traces
P.iters  = 4000;     % training iterations
P.bs     = 64;       % batch size
P.lr     = 1e-3;     % Adam learning rate (start value)
P.lrEnd  = 1e-5;     % final learning rate for the cosine schedule
P.lamEI  = 1.0;      % EI weight
P.gmax   = 80;       % max shift (samples), non-circular
P.tdeep  = 0.6;      % "deep window" for metrics (s)
P.tshallow = 0.3;    % "shallow window" for metrics (s)
P.hfBand = [40 70];  % band (Hz) for the descriptive high-frequency statistics
P.evalEvery = 500;   % iterations between test-set SNR records (monitoring only)

% Quick run for debugging: set environment variable EI_SMOKE=1 before starting MATLAB
if ~isempty(getenv('EI_SMOKE'))
    P.iters = 30; P.ntr = 200; P.nte = 50; P.evalEvery = 10;
end
P.seed  = envNum('EI_SEED', 0);
P.Qerr  = envNum('EI_QERR', P.Qerr);
P.iters = envNum('EI_ITERS', P.iters);
P.cfg   = envNum('EI_CFG', [1 2 3]);           % which of the three self-supervised variants to train
P.sched = envStr('EI_SCHED', 'cosine');
P.scen  = envStr('EI_SCEN', 'stationary');     % reflectivity statistics, see reflectivity()
P.frefNyq = 1 / (2 * P.dt);
P.fref  = envNum('EI_FREF', P.frefNyq);        % dispersion reference frequency (Hz)
P.nva  = 500;                                 % held-out traces for label-free validation losses
rng(P.seed);

t  = (0:P.N-1)' * P.dt;
Qm = P.Q * (1 + P.Qerr);                       % Q assumed by the method

A_true  = attenMatrix(P.Q, P.N, P.dt, P.fref); % generates data
A_model = attenMatrix(Qm,  P.N, P.dt, P.fref); % used inside the losses
useGPU  = canUseGPU;
if useGPU, A_model = gpuArray(single(A_model)); else, A_model = single(A_model); end

%% ---------------- synthetic data ----------------
[Xtr, Ytr] = makeData(P.ntr, A_true, P);
[Xte, Yte] = makeData(P.nte, A_true, P);
[~,   Yva] = makeData(P.nva, A_true, P);         % labels discarded: validation uses y only
s = std(Ytr(:));                               % input scaling for the network
deep = t > P.tdeep;
shallow = t < P.tshallow;

%% ---------------- baseline: Tikhonov (oracle-tuned mu) ----------------
best = -inf;
for mu = logspace(-6, 0, 25)
    Xh = (A_model' * A_model + mu * eye(P.N, 'like', A_model)) \ (A_model' * single(Yte));
    Xh = double(gather(Xh));
    v = snrdb(Xte, Xh);
    if v > best, best = v; Xtik = Xh; mubest = mu; end
end
fprintf('Tikhonov (mu=%.1e): SNR all %.2f dB | deep %.2f dB\n', mubest, ...
        snrdb(Xte, Xtik), snrdb(Xte(deep,:), Xtik(deep,:)));

%% ---------------- f_max(t) band-limit operator ----------------
Bf = fmaxMatrix(t, P.dt, Qm, P.snrdB, 40);
if useGPU, Bf = gpuArray(single(Bf)); else, Bf = single(Bf); end

%% ---------------- train three self-supervised variants ----------------
cfgs = struct('name', {'MC only', 'MC + EI', 'MC + EI + fmax'}, ...
              'lam',  {0, P.lamEI, P.lamEI}, ...
              'useB', {false, false, true});
R = struct();
for k = 1:numel(P.cfg)
    c = P.cfg(k);
    rng(1000 * P.seed + c);                    % each variant independent of run order and subset
    net = buildNet();
    [net, hist] = trainEI(net, Ytr, Xte, Yte, Yva, deep, s, t, A_model, Bf, cfgs(c), P, useGPU);
    Xh  = predictX(net, Yte, s, t, Bf, cfgs(c).useB, useGPU);
    R(k).name = cfgs(c).name;
    R(k).cfg  = c;
    R(k).snr_all  = snrdb(Xte, Xh);
    R(k).snr_deep = snrdb(Xte(deep,:), Xh(deep,:));
    R(k).Xh = Xh;
    R(k).hist = hist;   % see histNames below; monitoring only
    fprintf('%-16s: SNR all %.2f dB | deep %.2f dB\n', R(k).name, R(k).snr_all, R(k).snr_deep);
end

%% ---------------- plot deep-window spectra ----------------
nd = nnz(deep); fz = (0:floor(nd/2))' / (nd * P.dt);
spec = @(X) mean(abs(fftRows(X(deep,:), nd)), 2);
figure('Color','w'); hold on;
plot(fz, spec(Xte), 'k', 'LineWidth', 2);
plot(fz, spec(Yte), 'Color', [.75 .75 .75]);
plot(fz, spec(Xtik), 'Color', [.5 .5 .5]);
cols = [0.82 0.55 0.17; 0.17 0.42 0.82; 0.16 0.62 0.44];
for k = 1:numel(R), plot(fz, spec(R(k).Xh), 'Color', cols(R(k).cfg,:), 'LineWidth', 1.3); end
legend([{'True','Attenuated','Tikhonov'}, {R.name}], 'Box', 'off');
xlim([0 90]); xlabel('Frequency (Hz)'); ylabel('Mean amplitude');
title(sprintf('Deep window (t > %.1f s), Q=%g, Q_{model}=%g, SNR=%g dB', P.tdeep, P.Q, Qm, P.snrdB));

%% ---------------- save outputs ----------------
tag = sprintf('Qerr%+.2f_snr%g_seed%d_it%d_%s', P.Qerr, P.snrdB, P.seed, P.iters, P.sched);
if ~strcmp(P.scen, 'stationary'), tag = [tag '_' P.scen]; end
if P.fref ~= P.frefNyq, tag = [tag sprintf('_fref%g', P.fref)]; end
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(gcf, fullfile(outDir, ['spectra_' tag '.png']), 'Resolution', 200);

ns = nnz(shallow); fzs = (0:floor(ns/2))' / (ns * P.dt);
specS = @(X) mean(abs(fftRows(X(shallow,:), ns)), 2);
figure('Color','w'); hold on;
plot(fzs, specS(Xte), 'k', 'LineWidth', 2);
plot(fzs, specS(Yte), 'Color', [.75 .75 .75]);
plot(fzs, specS(Xtik), 'Color', [.5 .5 .5]);
for k = 1:numel(R), plot(fzs, specS(R(k).Xh), 'Color', cols(R(k).cfg,:), 'LineWidth', 1.3); end
legend([{'True','Attenuated','Tikhonov'}, {R.name}], 'Box', 'off');
xlim([0 90]); xlabel('Frequency (Hz)'); ylabel('Mean amplitude');
title(sprintf('Shallow window (t < %.1f s), scenario %s', P.tshallow, P.scen));
exportgraphics(gcf, fullfile(outDir, ['spectraShallow_' tag '.png']), 'Resolution', 200);
specMeans = struct('fDeep', fz, 'fShallow', fzs, 'trueDeep', spec(Xte), 'trueShallow', specS(Xte), ...
    'tikDeep', spec(Xtik), 'tikShallow', specS(Xtik), ...
    'estDeep', {cellfun(spec, {R.Xh}, 'UniformOutput', false)}, ...
    'estShallow', {cellfun(specS, {R.Xh}, 'UniformOutput', false)});

figure('Color','w'); hold on;
for k = 1:numel(R), plot(R(k).hist(:,1), R(k).hist(:,2), 'Color', cols(R(k).cfg,:), 'LineWidth', 1.3); end
yline(snrdb(Xte, Xtik), '--', 'Color', [.5 .5 .5]);
legend([{R.name}, {'Tikhonov (oracle \mu)'}], 'Box', 'off', 'Location', 'best');
xlabel('Iteration'); ylabel('Test SNR, all samples (dB)');
exportgraphics(gcf, fullfile(outDir, ['curves_' tag '.png']), 'Resolution', 200);

names = [{sprintf('Tikhonov(mu=%.1e)', mubest)}, {R.name}]';
snrAll  = [snrdb(Xte, Xtik);                 [R.snr_all]'];
snrDeep = [snrdb(Xte(deep,:), Xtik(deep,:)); [R.snr_deep]'];
allX = [{Xtik}, {R.Xh}];
snrShallow = cellfun(@(Z) snrdb(Xte(shallow,:), Z(shallow,:)), allX)';
[hfRatioDeep, snrHFDeep]       = cellfun(@(Z) bandStats(Xte, Z, deep, P), allX);
[hfRatioShallow, snrHFShallow] = cellfun(@(Z) bandStats(Xte, Z, shallow, P), allX);
hfRatioDeep = hfRatioDeep'; snrHFDeep = snrHFDeep'; hfRatioShallow = hfRatioShallow'; snrHFShallow = snrHFShallow';
writetable(table(names, snrAll, snrDeep, snrShallow, hfRatioDeep, snrHFDeep, hfRatioShallow, snrHFShallow), ...
    fullfile(outDir, ['summary_' tag '.csv']));
hists = {R.hist}; rname = {R.name};
save(fullfile(outDir, ['history_' tag '.mat']), 'P', 'hists', 'rname', 'mubest', 'specMeans');
histNames = {'iter','snrAll','snrDeep','loss','mc','ei','valMC','valEI','lr'};
H = cell(numel(R), 1);
for k = 1:numel(R)
    T = array2table(R(k).hist, 'VariableNames', histNames);
    T.config = repmat({R(k).name}, height(T), 1);
    H{k} = T;
end
writetable(vertcat(H{:}), fullfile(outDir, ['history_' tag '.csv']));
%% ======================================================================
%%                           local functions
%% ======================================================================
function A = attenMatrix(Q, N, dt, fh)
% Kolsky constant-Q operator (amplitude + dispersion). Column j = attenuated impulse at t_j.
% fh = dispersion reference frequency (Hz).
M  = 2 * N;                                   % zero-pad to avoid wrap-around
f  = (0:M/2)' / (M * dt); f(1) = 1e-6;
A  = zeros(N, N);
for j = 1:N
    tau = (j - 1) * dt;
    H = exp(-pi * f * tau / Q) .* exp(-1i * 2 * pi * f * tau .* (1 - log(f / fh) / (pi * Q)));
    Hfull = [H; conj(H(end-1:-1:2))];
    col = ifft(Hfull, 'symmetric');
    A(:, j) = col(1:N);
end
end

function [X, Y] = makeData(n, A, P)
tw = (-40:40)' * P.dt;
w  = (1 - 2 * (pi * P.fp * tw).^2) .* exp(-(pi * P.fp * tw).^2);
r  = reflectivity(n, P);
X  = conv2(r, w, 'same');                     % trace without attenuation
Y0 = A * X;
sig = sqrt(mean(Y0(:).^2)) * 10^(-P.snrdB / 20);
Y  = Y0 + sig * randn(size(Y0));
end

function r = reflectivity(n, P)
% Reflectivity series; u runs from 0 (first sample) to 1 (last sample).
% Non-stationary cases are rescaled to the stationary variance (0.08), so that
% loss magnitudes stay comparable. The stationary branch keeps its original draws.
N = P.N; u = (0:N-1)' / (N - 1);
switch P.scen
    case 'stationary'
        r = (rand(N, n) < 0.08) .* randn(N, n);
        return
    case 'density'       % dense shallow -> sparse deep, mean density 0.08
        r = (rand(N, n) < (0.14 - 0.12 * u)) .* randn(N, n);
    case 'amp'           % strong shallow -> weak deep
        r = (1.6 - 1.2 * u) .* (rand(N, n) < 0.08) .* randn(N, n);
    case {'broad2narrow', 'narrow2broad'}
        a = min(1, max(0, (u - 0.2) / 0.4));   % 0 above u = 0.2, 1 below u = 0.6
        if strcmp(P.scen, 'narrow2broad'), a = 1 - a; end
        rw = (rand(N, n) < 0.08) .* randn(N, n);   % white part
        rc = (rand(N, n) < 0.08) .* randn(N, n);   % part low-passed by a Gaussian, sigma 2 samples
        k  = (-6:6)'; g = exp(-k.^2 / (2 * 2^2));
        rc = conv2(rc, g, 'same');
        rc = rc * sqrt(0.08) / sqrt(mean(rc(:).^2));
        r  = sqrt(1 - a) .* rw + sqrt(a) .* rc;
    otherwise
        error('Unknown scenario %s', P.scen);
end
r = r * sqrt(0.08) / sqrt(mean(r(:).^2));
end

function [ratio, snrBand] = bandStats(X, Xh, win, P)
% Descriptive only: energy ratio and SNR inside P.hfBand within a time window.
n = nnz(win); f = (0:floor(n/2))' / (n * P.dt);
b = f >= P.hfBand(1) & f <= P.hfBand(2);
F  = fftRows(X(win,:), n);  Fh = fftRows(Xh(win,:), n);
ratio   = sum(abs(Fh(b,:)).^2, 'all') / sum(abs(F(b,:)).^2, 'all');
snrBand = 10 * log10(sum(abs(F(b,:)).^2, 'all') / sum(abs(F(b,:) - Fh(b,:)).^2, 'all'));
end

function B = fmaxMatrix(t, dt, Q, snrdB, L)
% Time-varying low-pass: f_max(t) = Q*ln(1/eps)/(pi*t), eps = noise/signal amplitude ratio.
% Row i = windowed-sinc kernel with cutoff f_max(t_i), applied around sample i.
N   = numel(t); fN = 1 / (2 * dt);
ep  = 10^(-snrdB / 20);
fmx = min(fN, max(10, Q * log(1 / ep) ./ (pi * max(t, dt))));
k   = (-L:L);
win = 0.5 * (1 + cos(pi * k / (L + 1)));
B   = zeros(N, N);
for i = 1:N
    fc = fmx(i) * dt;                          % normalized cutoff (cycles/sample)
    h  = 2 * fc * mysinc(2 * fc * k) .* win;
    h  = h / sum(h);
    idx = i + k; ok = idx >= 1 & idx <= N;
    B(i, idx(ok)) = h(ok);
end
end

function y = mysinc(x)
y = ones(size(x)); nz = x ~= 0;
y(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
end

function net = buildNet()
% Dilated 1D CNN; input channels: [scaled y ; normalized time]. Output: stationary trace.
C = 48; layers = [sequenceInputLayer(2, 'MinLength', 256, 'Normalization', 'none')];
for d = [1 2 4 8 16 1 2 4 8 16]
    layers = [layers
        convolution1dLayer(9, C, 'Padding', 'same', 'DilationFactor', d)
        leakyReluLayer(0.1)]; %#ok<AGROW>
end
layers = [layers; convolution1dLayer(1, 1)];
net = dlnetwork(layers);
end

function D = toCBT(M)          % N x B matrix -> 1 x B x N formatted dlarray
D = dlarray(permute(M, [3 2 1]), 'CBT');
end

function M = fromCBT(D)        % 1 x B x N dlarray -> N x B
M = permute(stripdims(D), [3 2 1]);
end

function Z = shiftNC(X, g)     % non-circular shift down by g samples (zero fill)
[n, b] = size(X);
if g >= 0, Z = [zeros(g, b, 'like', X); X(1:n-g, :)];
else,      Z = [X(1-g:n, :); zeros(-g, b, 'like', X)]; end
end

function X = applyNet(net, Ymat, s, Tch, Bf, useB)
X = fromCBT(forward(net, cat(1, toCBT(Ymat / s), Tch)));
if useB, X = Bf * X; end
end

function [loss, grads, lmc, lei] = modelLoss(net, Y, s, Tch, A, Bf, cfg, gmax)
Xh  = applyNet(net, Y, s, Tch, Bf, cfg.useB);
lmc = mean(sum((A * Xh - Y).^2, 1));
lei = dlarray(0);
if cfg.lam > 0
    g   = randi([-gmax gmax]);
    Z   = shiftNC(Xh, g);                     % T_g x^
    Xz  = applyNet(net, A * Z, s, Tch, Bf, cfg.useB);   % f(A T_g x^)
    lei = mean(sum((Xz - Z).^2, 1));
end
loss  = lmc + cfg.lam * lei;
grads = dlgradient(loss, net.Learnables);
end

function [net, hist] = trainEI(net, Ytr, Xte, Yte, Yva, deep, s, t, A, Bf, cfg, P, useGPU)
avg = []; avgSq = []; hist = [];
Tch = repmat(reshape(single(t / t(end)), 1, 1, []), 1, P.bs, 1);
Tch = dlarray(Tch, 'CBT'); if useGPU, Tch = gpuArray(Tch); end
t0 = tic;
for it = 1:P.iters
    idx = randi(size(Ytr, 2), 1, P.bs);
    Y = dlarray(single(Ytr(:, idx)));
    if useGPU, Y = gpuArray(Y); end
    [loss, grads, lmc, lei] = dlfeval(@modelLoss, net, Y, s, Tch, A, Bf, cfg, P.gmax);
    lr = lrAt(P, it);
    [net, avg, avgSq] = adamupdate(net, grads, avg, avgSq, it, lr);
    if mod(it, P.evalEvery) == 0 || it == P.iters
        % Test and validation data are only monitored: nothing below feeds back into training.
        Xh = predictX(net, Yte, s, t, Bf, cfg.useB, useGPU);
        row = [it, snrdb(Xte, Xh), snrdb(Xte(deep,:), Xh(deep,:)), ...
               double(extractdata(gather(loss))), double(extractdata(gather(lmc))), double(extractdata(gather(lei)))];
        [vm, ve] = valLosses(net, Yva, s, t, A, Bf, cfg.useB, useGPU);
        row = [row, vm, ve, lr]; %#ok<AGROW>
        hist = [hist; row]; %#ok<AGROW>
        fprintf('  [%s] it %d  loss %.4g  MC %.4g  EI %.4g | test SNR %.2f / deep %.2f dB (%.0f s)\n', cfg.name, it, ...
            row(4), row(5), row(6), row(2), row(3), toc(t0));
    end
end
end

function lr = lrAt(P, it)
if strcmp(P.sched, 'cosine')
    lr = P.lrEnd + 0.5 * (P.lr - P.lrEnd) * (1 + cos(pi * (it - 1) / max(1, P.iters - 1)));
else
    lr = P.lr;
end
end

function [vm, ve] = valLosses(net, Yva, s, t, A, Bf, useB, useGPU)
% Label-free validation: held-out MC term and EI term with fixed shifts.
Ad = double(gather(A));
Xh = predictX(net, Yva, s, t, Bf, useB, useGPU);
vm = mean(sum((Ad * Xh - Yva).^2, 1));
shifts = [-60 -20 20 60];
ve = 0;
for g = shifts
    Z  = shiftNC(Xh, g);
    Xz = predictX(net, Ad * Z, s, t, Bf, useB, useGPU);
    ve = ve + mean(sum((Xz - Z).^2, 1)) / numel(shifts);
end
end

function v = envNum(name, default)
v = default;
e = getenv(name);
if ~isempty(e), v = str2num(e); end %#ok<ST2NM>  % allows lists such as "1 2"
end

function v = envStr(name, default)
v = default;
e = getenv(name);
if ~isempty(e), v = e; end
end

function X = predictX(net, Y, s, t, Bf, useB, useGPU)
X = zeros(size(Y)); bs = 100;
for i0 = 1:bs:size(Y, 2)
    idx = i0:min(i0 + bs - 1, size(Y, 2));
    Tch = dlarray(repmat(reshape(single(t / t(end)), 1, 1, []), 1, numel(idx), 1), 'CBT');
    Yb  = dlarray(single(Y(:, idx)));
    if useGPU, Tch = gpuArray(Tch); Yb = gpuArray(Yb); end
    Xb  = fromCBT(predict(net, cat(1, toCBT(Yb / s), Tch)));
    if useB, Xb = Bf * Xb; end
    X(:, idx) = double(gather(extractdata(Xb)));
end
end

function v = snrdb(X, Xh)
v = 10 * log10(sum(X(:).^2) / sum((X(:) - Xh(:)).^2));
end

function S = fftRows(X, n)
F = fft(X, n, 1); S = F(1:floor(n/2) + 1, :);
end
