%% aggregate_nonstat.m -- applies the rules in PLAN_NONSTAT_ZH.md.
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
seeds = [1 2 3]; iters = 4000;
if ~isempty(getenv('EI_AGG_ITERS')), iters = str2double(getenv('EI_AGG_ITERS')); end  % smoke test only
if ~isempty(getenv('EI_AGG_SEEDS')), seeds = str2num(getenv('EI_AGG_SEEDS')); end     %#ok<ST2NM> % smoke test only
scens = {'density', 'amp', 'broad2narrow', 'narrow2broad'};
names = {'Tikhonov', 'MC only', 'MC + EI'};
extra = {'snrShallow', 'hfRatioDeep', 'snrHFDeep', 'hfRatioShallow', 'snrHFShallow'};

rows = {};
for sc = [{'stationary'}, scens]
    for s = seeds
        if strcmp(sc{1}, 'stationary')   % reused baseline, always the real 4000-step runs
            base = sprintf('Qerr+0.00_snr30_seed%d_it4000_cosine', s);
        else
            base = sprintf('Qerr+0.00_snr30_seed%d_it%d_cosine_%s', s, iters, sc{1});
        end
        T = readtable(fullfile(outDir, ['summary_' base '.csv']), 'TextType', 'char');
        H = readtable(fullfile(outDir, ['history_' base '.csv']), 'TextType', 'char');
        for n = 1:numel(names)
            if n == 1, r = T(startsWith(T.names, 'Tikhonov'), :);
            else,      r = T(strcmp(T.names, names{n}), :); end
            ex = nan(1, numel(extra));
            for e = 1:numel(extra)
                if ismember(extra{e}, T.Properties.VariableNames), ex(e) = r.(extra{e}); end
            end
            best = nan;
            if n > 1, best = max(H.snrAll(strcmp(H.config, names{n}))); end
            rows(end+1, :) = [{sc{1}, s, names{n}, r.snrAll, r.snrDeep, best}, num2cell(ex)]; %#ok<SAGROW>
        end
    end
end
R = cell2table(rows, 'VariableNames', [{'scenario','seed','method','snrAll','snrDeep','bestAll'}, extra]);
writetable(R, fullfile(outDir, 'nonstat_table.csv'));
disp(R);

get = @(sc, m, col) arrayfun(@(s) R.(col)(strcmp(R.scenario, sc) & R.seed == s & strcmp(R.method, m)), seeds);
fid = fopen(fullfile(outDir, 'nonstat_verdicts.txt'), 'w');
pr = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);

d0 = get('stationary', 'MC + EI', 'snrAll') - get('stationary', 'MC only', 'snrAll');
pr('Stationary baseline (reused): MC+EI - MC, all = %s (mean %.2f); deep = %s (mean %.2f); best recorded = %s (mean %.2f)\n', ...
    mat2str(d0, 3), mean(d0), ...
    mat2str(get('stationary', 'MC + EI', 'snrDeep') - get('stationary', 'MC only', 'snrDeep'), 3), ...
    mean(get('stationary', 'MC + EI', 'snrDeep') - get('stationary', 'MC only', 'snrDeep')), ...
    mat2str(get('stationary', 'MC + EI', 'bestAll') - get('stationary', 'MC only', 'bestAll'), 3), ...
    mean(get('stationary', 'MC + EI', 'bestAll') - get('stationary', 'MC only', 'bestAll')));

e1 = cell(1, numel(scens));
for c = 1:numel(scens)
    sc = scens{c};
    pr('\n--- %s ---\n', sc);
    d = get(sc, 'MC + EI', 'snrAll') - get(sc, 'MC only', 'snrAll'); e1{c} = rule(d);
    pr('E1 MC+EI - MC, all  = %s (mean %.2f) -> %s\n', mat2str(d, 3), mean(d), e1{c});
    d = get(sc, 'MC + EI', 'snrDeep') - get(sc, 'MC only', 'snrDeep');
    pr('E2 MC+EI - MC, deep = %s (mean %.2f) -> %s\n', mat2str(d, 3), mean(d), rule(d));
    d = get(sc, 'MC + EI', 'snrAll') - get(sc, 'Tikhonov', 'snrAll');
    pr('E3 MC+EI - Tikhonov, all = %s (mean %.2f) -> %s\n', mat2str(d, 3), mean(d), rule(d));
    d = get(sc, 'MC + EI', 'bestAll') - get(sc, 'MC only', 'bestAll');
    pr('   descriptive: best recorded MC+EI - MC = %s (mean %.2f)\n', mat2str(d, 3), mean(d));
    for m = names
        pr('   descriptive %-9s: all %.2f, deep %.2f, shallow %.2f | HF ratio deep %.2f, shallow %.2f | HF SNR deep %.2f, shallow %.2f\n', m{1}, ...
            mean(get(sc, m{1}, 'snrAll')), mean(get(sc, m{1}, 'snrDeep')), mean(get(sc, m{1}, 'snrShallow')), ...
            mean(get(sc, m{1}, 'hfRatioDeep')), mean(get(sc, m{1}, 'hfRatioShallow')), ...
            mean(get(sc, m{1}, 'snrHFDeep')), mean(get(sc, m{1}, 'snrHFShallow')));
    end
end

kept = strcmp(e1, 'advantage retained');
spectral = ismember(scens, {'broad2narrow', 'narrow2broad'});
if all(kept), overall = 'EI premise robust to these statistical changes';
elseif ~any(kept(spectral)), overall = 'EI advantage depends on spectral stationarity';
else, overall = 'partly sensitive, see scenarios'; end
pr('\nOverall (E1): %s\n', overall);
fclose(fid);

function v = rule(d)
if all(d > 0) && mean(d) >= 1, v = 'advantage retained';
elseif all(d > 0), v = 'positive but below 1 dB';
elseif all(d < 0) && mean(d) <= -1, v = 'harmful';
else, v = 'undecided'; end
end
