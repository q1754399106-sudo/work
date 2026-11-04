%% decompose_qerr.m -- post hoc, descriptive. Splits the Tikhonov loss under a Q error into
% an amplitude part and a dispersion part, by giving the operator a wrong Q in only one term.
% Uses the same test traces as ei_qcomp_net.m for seeds 1-3 (same rng order). No training.
clear;
P.N = 256; P.dt = 0.004; P.Q = 50; P.fp = 25; P.snrdB = 30; P.ntr = 4000; P.nte = 500;
t = (0:P.N-1)' * P.dt; deep = t > 0.6;
Atrue = atten2(P.Q, P.Q, P.N, P.dt);
cases = {'Q correct', 50, 50; ...
         'Q=60 both', 60, 60; 'Q=60 amplitude only', 60, 50; 'Q=60 dispersion only', 50, 60; ...
         'Q=40 both', 40, 40; 'Q=40 amplitude only', 40, 50; 'Q=40 dispersion only', 50, 40};
seeds = [1 2 3];
out = zeros(size(cases, 1), numel(seeds), 2);
for is = 1:numel(seeds)
    rng(seeds(is));
    makeData(P.ntr, Atrue, P);                 % consume the training draw, as the main script does
    [Xte, Yte] = makeData(P.nte, Atrue, P);
    for c = 1:size(cases, 1)
        A = atten2(cases{c, 2}, cases{c, 3}, P.N, P.dt);
        best = -inf;
        for mu = logspace(-6, 0, 25)
            Xh = (A' * A + mu * eye(P.N)) \ (A' * Yte);
            v = snrdb(Xte, Xh);
            if v > best, best = v; out(c, is, 1) = v; out(c, is, 2) = snrdb(Xte(deep, :), Xh(deep, :)); end
        end
    end
end
names = cases(:, 1);
T = table(names, mean(out(:, :, 1), 2), mean(out(:, :, 2), 2), 'VariableNames', {'case', 'snrAll', 'snrDeep'});
T.lossAll = T.snrAll(1) - T.snrAll;
T.lossDeep = T.snrDeep(1) - T.snrDeep;
disp(T);
writetable(T, fullfile(fileparts(mfilename('fullpath')), 'outputs', 'qerr_decomposition.csv'));

function A = atten2(Qamp, Qphase, N, dt)
M  = 2 * N;
f  = (0:M/2)' / (M * dt); f(1) = 1e-6;
fh = 1 / (2 * dt);
A  = zeros(N, N);
for j = 1:N
    tau = (j - 1) * dt;
    H = exp(-pi * f * tau / Qamp) .* exp(-1i * 2 * pi * f * tau .* (1 - log(f / fh) / (pi * Qphase)));
    Hfull = [H; conj(H(end-1:-1:2))];
    col = ifft(Hfull, 'symmetric');
    A(:, j) = col(1:N);
end
end

function [X, Y] = makeData(n, A, P)
tw = (-40:40)' * P.dt;
w  = (1 - 2 * (pi * P.fp * tw).^2) .* exp(-(pi * P.fp * tw).^2);
r  = (rand(P.N, n) < 0.08) .* randn(P.N, n);
X  = conv2(r, w, 'same');
Y0 = A * X;
sig = sqrt(mean(Y0(:).^2)) * 10^(-P.snrdB / 20);
Y  = Y0 + sig * randn(size(Y0));
end

function v = snrdb(X, Xh)
v = 10 * log10(sum(X(:).^2) / sum((X(:) - Xh(:)).^2));
end
