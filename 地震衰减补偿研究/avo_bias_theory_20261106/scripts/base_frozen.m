function base_frozen(exdir, srcdir, inputsdir, pattern)
% Frozen preset, offset domain: prepare_models_v2, uq_select, unchanged; uq_build_prior repeated below with the known sigma instead of the WindowMAD estimate
% (the archive's KnownNoise mode; in dense scenes the 0.02-0.25 s window holds signal, and all other methods here are given sigma). The Q-estimation step is skipped: Q is given, as for all methods.
addpath(srcdir); maxNumCompThreads(4); txt = fileread(fullfile(inputsdir, 'uq_spec.json')); o = jsondecode(txt); ix = find(strcmp({o.presets.name}, 'Frozen'), 1); preset = o.presets(ix); fs = dir(fullfile(exdir, [pattern '_meta.json'])); c0 = config_model();
for k = 1:numel(fs)
 m = jsondecode(fileread(fullfile(exdir, fs(k).name))); out = fullfile(exdir, [m.name '_FrozenSplit.bin']); if isfile(out), continue; end
 D = rd(fullfile(exdir, [m.name '_D.bin'])); D = reshape(D, 601, 11, m.nc); c = c0; c.Q = m.Q(:)'; c.fpeak = m.fp; t0 = tic; models = prepare_models_v2(m.offsets(:)', c); t1 = toc(t0); S = zeros(size(D)); Dd = S; ep = 0;
 for q = 1:m.nc
  prior = known_prior(D(:, :, q), models, c, preset, m.sig(:, q)'); r = uq_select(D(:, :, q), models, prior, preset); S(:, :, q) = r.split; Dd(:, :, q) = r.direct; ep = ep + r.endpoint;
 end
 wr(out, S); wr(fullfile(exdir, [m.name '_FrozenDirect.bin']), Dd); fprintf('%s Frozen models %.1fs total %.1fs endpoint %d\n', m.name, t1, toc(t0), ep);
end
end
function prior = known_prior(D, models, c, preset, sigma)
% uq_build_prior with domain 'Offset' and the given sigma; nothing else changed.
[n, nx] = size(D); lambda = zeros(1, nx); pilot = zeros(n, nx);
for j = 1:nx
 [~, lambda(j)] = discrepancy_apply(models{j}.direct, D(:, j), c.lambda_grid, sqrt(n) * sigma(j));
 F = models{j}.direct; u = F.U' * D(:, j); pilot(:, j) = F.V * (F.z ./ (F.z.^2 + preset.pilot_multiplier * lambda(j)) .* u);
end
H = sparse(n * nx, n * nx);
if preset.gamma > 0
 [H, ~] = learn_structure_operator(pilot, c, structure_options()); R = avo_projector_v2(models); H = H * R;
end
prior = struct('H', H, 'sigma', sigma, 'lambda', lambda, 'pilot', pilot);
end
function x = rd(p), f = fopen(p, 'r'); x = fread(f, inf, 'double'); fclose(f); end
function wr(p, x), f = fopen(p, 'w'); fwrite(f, x, 'double'); fclose(f); end
