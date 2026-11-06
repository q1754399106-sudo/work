function base_fcr(exdir, srcdir, pattern)
% FullCalRefit: reflectivity_refit_prepared, unchanged, one dictionary per offset. Known sigma, given Q, dictionary at the given frequency.
addpath(srcdir); maxNumCompThreads(4); fs = dir(fullfile(exdir, [pattern '_meta.json'])); cm0 = config_model();
for k = 1:numel(fs)
 m = jsondecode(fileread(fullfile(exdir, fs(k).name))); out = fullfile(exdir, [m.name '_FCR.bin']); if isfile(out), continue; end
 D = rd(fullfile(exdir, [m.name '_D.bin'])); D = reshape(D, 601, 11, m.nc); S = zeros(size(D)); cm = cm0; cm.Q = m.Q(:)'; cm.fpeak = m.fp; t0 = tic; flags = 0;
 for j = 1:11
  d = prune_dictionary(cm, m.offsets(j)); dicts = repmat({d}, 1, m.nc); sg = reshape(m.sig(j, :), 1, []); Dj = reshape(D(:, j, :), 601, m.nc);
  r = reflectivity_refit_prepared(Dj, dicts, sg); S(:, j, :) = reshape(r.gather, 601, 1, m.nc); flags = flags + r.optimization_flags;
 end
 wr(out, S); fprintf('%s FCR %.1fs flags %d\n', m.name, toc(t0), flags);
end
end
function x = rd(p), f = fopen(p, 'r'); x = fread(f, inf, 'double'); fclose(f); end
function wr(p, x), f = fopen(p, 'w'); fwrite(f, x, 'double'); fclose(f); end
