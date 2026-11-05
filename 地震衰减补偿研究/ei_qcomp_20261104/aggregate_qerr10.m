%% aggregate_qerr10.m -- applies the rules in PLAN_QERR10_ZH.md and lists the five Q values.
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
seeds = [1 2 3]; iters = 4000;
if ~isempty(getenv('EI_AGG_ITERS')), iters = str2double(getenv('EI_AGG_ITERS')); end  % smoke test only
if ~isempty(getenv('EI_AGG_SEEDS')), seeds = str2num(getenv('EI_AGG_SEEDS')); end     %#ok<ST2NM> % smoke test only
qs = {'-0.20', '-0.10', '+0.00', '+0.10', '+0.20'};
fresh = {'-0.10', '+0.10'};
names = {'Tikhonov', 'MC only', 'MC + EI'};

rows = {};
for q = 1:numel(qs)
    it = 4000; if ismember(qs{q}, fresh), it = iters; end   % reused runs are always the real 4000-step ones
    Qm = 50 * (1 + str2double(qs{q}));
    for s = seeds
        base = sprintf('Qerr%s_snr30_seed%d_it%d_cosine', qs{q}, s, it);
        T = readtable(fullfile(outDir, ['summary_' base '.csv']), 'TextType', 'char');
        H = readtable(fullfile(outDir, ['history_' base '.csv']), 'TextType', 'char');
        for n = 1:numel(names)
            if n == 1, r = T(startsWith(T.names, 'Tikhonov'), :); best = nan; vm = nan; ve = nan;
            else
                r = T(strcmp(T.names, names{n}), :);
                h = H(strcmp(H.config, names{n}), :);
                best = max(h.snrAll); vm = h.valMC(end); ve = h.valEI(end);
            end
            rows(end+1, :) = {qs{q}, Qm, s, names{n}, r.snrAll, r.snrDeep, best, vm, ve}; %#ok<SAGROW>
        end
    end
end
R = cell2table(rows, 'VariableNames', {'Qerr','Qmodel','seed','method','snrAll','snrDeep','bestAll','valMC','valEI'});
writetable(R, fullfile(outDir, 'qerr10_table.csv'));
disp(R);

get = @(q, m, col) arrayfun(@(s) R.(col)(strcmp(R.Qerr, q) & R.seed == s & strcmp(R.method, m)), seeds);
fid = fopen(fullfile(outDir, 'qerr10_verdicts.txt'), 'w');
pr = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);

verd = {};
for q = fresh
    d = get(q{1}, 'MC + EI', 'snrAll') - get(q{1}, 'MC only', 'snrAll'); verd{end+1} = rule(d); %#ok<SAGROW>
    pr('F1 Qerr %s: MC+EI - MC = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), verd{end});
    d = get(q{1}, 'MC + EI', 'snrAll') - get(q{1}, 'Tikhonov', 'snrAll'); verd{end+1} = rule(d); %#ok<SAGROW>
    pr('F2 Qerr %s: MC+EI - Tikhonov = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), verd{end});
end
if all(strcmp(verd, 'advantage retained')), overall = 'EI keeps its advantage with Q within 10 percent';
else, overall = 'EI needs Q better than 10 percent, see the lines above'; end
pr('Overall: %s\n\nDescriptive, mean over seeds\n', overall);

pr('%-16s %s\n', 'Q model', sprintf('%8g', 50 * (1 + str2double(qs))));
for m = names
    all_ = cellfun(@(q) mean(get(q, m{1}, 'snrAll')), qs);
    deep = cellfun(@(q) mean(get(q, m{1}, 'snrDeep')), qs);
    ref  = all_(strcmp(qs, '+0.00'));
    pr('%-16s %s  (all)\n',  m{1}, sprintf('%8.2f', all_));
    pr('%-16s %s  (deep)\n', '',   sprintf('%8.2f', deep));
    pr('%-16s %s  (loss vs Q correct, all)\n', '', sprintf('%8.2f', ref - all_));
end
for m = names(2:3)
    pr('%-16s %s  (best recorded, all)\n', m{1}, sprintf('%8.2f', cellfun(@(q) mean(get(q, m{1}, 'bestAll')), qs)));
    pr('%-16s %s  (held-out MC)\n', '', sprintf('%8.4f', cellfun(@(q) mean(get(q, m{1}, 'valMC')), qs)));
    pr('%-16s %s  (held-out EI)\n', '', sprintf('%8.4f', cellfun(@(q) mean(get(q, m{1}, 'valEI')), qs)));
end
pr('%-16s %s  (MC+EI - MC, all)\n', '', sprintf('%8.2f', cellfun(@(q) mean(get(q, 'MC + EI', 'snrAll') - get(q, 'MC only', 'snrAll')), qs)));
pr('%-16s %s  (MC+EI - Tikhonov, all)\n', '', sprintf('%8.2f', cellfun(@(q) mean(get(q, 'MC + EI', 'snrAll') - get(q, 'Tikhonov', 'snrAll')), qs)));
fclose(fid);

function v = rule(d)
if all(d > 0) && mean(d) >= 1, v = 'advantage retained';
elseif all(d > 0), v = 'positive but below 1 dB';
elseif all(d < 0) && mean(d) <= -1, v = 'harmful';
else, v = 'undecided'; end
end
