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
% Python/numpy linear prototype: ei_qcomp_linear.py

clear; rng(0);

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
P.lr     = 1e-3;     % Adam learning rate
P.lamEI  = 1.0;      % EI weight
P.gmax   = 80;       % max shift (samples), non-circular
P.tdeep  = 0.6;      % "deep window" for metrics (s)
P.evalEvery = 500;   % iterations between test-set SNR records (monitoring only)

% Quick run for debugging: set environment variable EI_SMOKE=1 before starting MATLAB
if ~isempty(getenv('EI_SMOKE'))
    P.iters = 30; P.ntr = 200; P.nte = 50; P.evalEvery = 10;
end

t  = (0:P.N-1)' * P.dt;
Qm = P.Q * (1 + P.Qerr);                       % Q assumed by the method

A_true  = attenMatrix(P.Q, P.N, P.dt);         % generates data
A_model = attenMatrix(Qm,  P.N, P.dt);         % used inside the losses
useGPU  = canUseGPU;
if useGPU, A_model = gpuArray(single(A_model)); else, A_model = single(A_model); end

%% ---------------- synthetic data ----------------
[Xtr, Ytr] = makeData(P.ntr, A_true, P);
[Xte, Yte] = makeData(P.nte, A_true, P);
s = std(Ytr(:));                               % input scaling for the network
deep = t > P.tdeep;

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
for c = 1:numel(cfgs)
    net = buildNet();
    [net, hist] = trainEI(net, Ytr, Xte, Yte, deep, s, t, A_model, Bf, cfgs(c), P, useGPU);
    Xh  = predictX(net, Yte, s, t, Bf, cfgs(c).useB, useGPU);
    R(c).name = cfgs(c).name;
    R(c).snr_all  = snrdb(Xte, Xh);
    R(c).snr_deep = snrdb(Xte(deep,:), Xh(deep,:));
    R(c).Xh = Xh;
    R(c).hist = hist;   % [iter, test SNR all, test SNR deep, loss, MC, EI]; monitoring only
    fprintf('%-16s: SNR all %.2f dB | deep %.2f dB\n', R(c).name, R(c).snr_all, R(c).snr_deep);
end

%% ---------------- plot deep-window spectra ----------------
nd = nnz(deep); fz = (0:floor(nd/2))' / (nd * P.dt);
spec = @(X) mean(abs(fftRows(X(deep,:), nd)), 2);
figure('Color','w'); hold on;
plot(fz, spec(Xte), 'k', 'LineWidth', 2);
plot(fz, spec(Yte), 'Color', [.75 .75 .75]);
plot(fz, spec(Xtik), 'Color', [.5 .5 .5]);
cols = [0.82 0.55 0.17; 0.17 0.42 0.82; 0.16 0.62 0.44];
for c = 1:numel(R), plot(fz, spec(R(c).Xh), 'Color', cols(c,:), 'LineWidth', 1.3); end
legend([{'True','Attenuated','Tikhonov'}, {R.name}], 'Box', 'off');
xlim([0 90]); xlabel('Frequency (Hz)'); ylabel('Mean amplitude');
title(sprintf('Deep window (t > %.1f s), Q=%g, Q_{model}=%g, SNR=%g dB', P.tdeep, P.Q, Qm, P.snrdB));

%% ---------------- save outputs ----------------
tag = sprintf('Qerr%+.2f_snr%g_seed0_it%d', P.Qerr, P.snrdB, P.iters);
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if ~isfolder(outDir), mkdir(outDir); end
exportgraphics(gcf, fullfile(outDir, ['spectra_' tag '.png']), 'Resolution', 200);

figure('Color','w'); hold on;
for c = 1:numel(R), plot(R(c).hist(:,1), R(c).hist(:,2), 'Color', cols(c,:), 'LineWidth', 1.3); end
yline(snrdb(Xte, Xtik), '--', 'Color', [.5 .5 .5]);
legend([{R.name}, {'Tikhonov (oracle \mu)'}], 'Box', 'off', 'Location', 'best');
xlabel('Iteration'); ylabel('Test SNR, all samples (dB)');
exportgraphics(gcf, fullfile(outDir, ['curves_' tag '.png']), 'Resolution', 200);

names = [{sprintf('Tikhonov(mu=%.1e)', mubest)}, {R.name}]';
snrAll  = [snrdb(Xte, Xtik);                 [R.snr_all]'];
snrDeep = [snrdb(Xte(deep,:), Xtik(deep,:)); [R.snr_deep]'];
writetable(table(names, snrAll, snrDeep), fullfile(outDir, ['summary_' tag '.csv']));
hists = {R.hist}; rname = {R.name};
save(fullfile(outDir, ['history_' tag '.mat']), 'P', 'hists', 'rname', 'mubest');
%% ======================================================================
%%                           local functions
%% ======================================================================
function A = attenMatrix(Q, N, dt)
% Kolsky constant-Q operator (amplitude + dispersion). Column j = attenuated impulse at t_j.
M  = 2 * N;                                   % zero-pad to avoid wrap-around
f  = (0:M/2)' / (M * dt); f(1) = 1e-6;
fh = 1 / (2 * dt);                            % reference frequency = Nyquist
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
r  = (rand(P.N, n) < 0.08) .* randn(P.N, n);
X  = conv2(r, w, 'same');                     % stationary trace
Y0 = A * X;
sig = sqrt(mean(Y0(:).^2)) * 10^(-P.snrdB / 20);
Y  = Y0 + sig * randn(size(Y0));
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

function [net, hist] = trainEI(net, Ytr, Xte, Yte, deep, s, t, A, Bf, cfg, P, useGPU)
avg = []; avgSq = []; hist = [];
Tch = repmat(reshape(single(t / t(end)), 1, 1, []), 1, P.bs, 1);
Tch = dlarray(Tch, 'CBT'); if useGPU, Tch = gpuArray(Tch); end
t0 = tic;
for it = 1:P.iters
    idx = randi(size(Ytr, 2), 1, P.bs);
    Y = dlarray(single(Ytr(:, idx)));
    if useGPU, Y = gpuArray(Y); end
    [loss, grads, lmc, lei] = dlfeval(@modelLoss, net, Y, s, Tch, A, Bf, cfg, P.gmax);
    [net, avg, avgSq] = adamupdate(net, grads, avg, avgSq, it, P.lr);
    if mod(it, P.evalEvery) == 0 || it == P.iters
        % Test data are only monitored here: nothing below feeds back into training.
        Xh = predictX(net, Yte, s, t, Bf, cfg.useB, useGPU);
        row = [it, snrdb(Xte, Xh), snrdb(Xte(deep,:), Xh(deep,:)), ...
               double(extractdata(gather(loss))), double(extractdata(gather(lmc))), double(extractdata(gather(lei)))];
        hist = [hist; row]; %#ok<AGROW>
        fprintf('  [%s] it %d  loss %.4g  MC %.4g  EI %.4g | test SNR %.2f / deep %.2f dB (%.0f s)\n', cfg.name, it, ...
            row(4), row(5), row(6), row(2), row(3), toc(t0));
    end
end
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
