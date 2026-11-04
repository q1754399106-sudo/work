%% aggregate_qerr.m -- applies the rules in PLAN_QERR_ZH.md.
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
seeds = [1 2 3]; qs = {'+0.00', '+0.20', '-0.20'}; iters = 4000;
if ~isempty(getenv('EI_AGG_ITERS')), iters = str2double(getenv('EI_AGG_ITERS')); end  % smoke test only
names = {'Tikhonov', 'MC only', 'MC + EI', 'MC + EI + fmax'};
rows = {};
for q = 1:numel(qs)
    for s = seeds
        T = readtable(fullfile(outDir, sprintf('summary_Qerr%s_snr30_seed%d_it%d_cosine.csv', qs{q}, s, iters)), 'TextType', 'char');
        for n = 1:numel(names)
            r = T(startsWith(T.names, names{n}) & (strcmp(names{n}, 'Tikhonov') | strcmp(T.names, names{n})), :);
            rows(end+1, :) = {qs{q}, s, names{n}, r.snrAll, r.snrDeep}; %#ok<SAGROW>
        end
    end
end
R = cell2table(rows, 'VariableNames', {'Qerr','seed','method','snrAll','snrDeep'});
writetable(R, fullfile(outDir, 'qerr_table.csv'));
disp(unstack(R(:, {'Qerr','seed','method','snrAll'}), 'snrAll', 'method'));
disp(unstack(R(:, {'Qerr','seed','method','snrDeep'}), 'snrDeep', 'method'));
get = @(q, m, col) arrayfun(@(s) R.(col)(strcmp(R.Qerr, q) & R.seed == s & strcmp(R.method, m)), seeds);
fid = fopen(fullfile(outDir, 'qerr_verdicts.txt'), 'w');
pr = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
for q = {'+0.20', '-0.20'}
    d = get(q{1}, 'MC + EI', 'snrAll') - get(q{1}, 'MC only', 'snrAll');
    pr('D1 Qerr %s: MC+EI - MC = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), rule(d));
    d = get(q{1}, 'MC + EI', 'snrAll') - get(q{1}, 'Tikhonov', 'snrAll');
    pr('D2 Qerr %s: MC+EI - Tikhonov = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), rule(d));
    d = get(q{1}, 'MC + EI + fmax', 'snrAll') - get(q{1}, 'MC + EI', 'snrAll');
    if all(d > 0) && mean(d) >= 0.5, v = 'gain';
    elseif all(d < 0) && mean(d) <= -0.5, v = 'harmful';
    else, v = 'no gain seen'; end
    pr('D3 Qerr %s: fmax - MC+EI = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), v);
    for m = names
        l = get('+0.00', m{1}, 'snrAll') - get(q{1}, m{1}, 'snrAll');
        ld = get('+0.00', m{1}, 'snrDeep') - get(q{1}, m{1}, 'snrDeep');
        pr('   loss vs Qerr 0, %s, %s: all %s (mean %.2f), deep %s (mean %.2f)\n', q{1}, m{1}, mat2str(l, 3), mean(l), mat2str(ld, 3), mean(ld));
    end
end
fclose(fid);

function v = rule(d)
if all(d > 0) && mean(d) >= 1, v = 'advantage retained';
elseif all(d > 0), v = 'positive but below 1 dB';
elseif all(d < 0) && mean(d) <= -1, v = 'harmful';
else, v = 'undecided'; end
end
