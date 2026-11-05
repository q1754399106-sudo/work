%% aggregate_fref.m -- applies the rules in PLAN_FREF_ZH.md.
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
seeds = [1 2 3]; iters = 4000;
if ~isempty(getenv('EI_AGG_ITERS')), iters = str2double(getenv('EI_AGG_ITERS')); end  % smoke test only
if ~isempty(getenv('EI_AGG_SEEDS')), seeds = str2num(getenv('EI_AGG_SEEDS')); end     %#ok<ST2NM> % smoke test only
qs = {'-0.10', '+0.00', '+0.10'};
errs = {'-0.10', '+0.10'};
refs = {'Nyquist', '25Hz'};
names = {'Tikhonov', 'MC only', 'MC + EI'};

rows = {};
for f = 1:numel(refs)
    for q = 1:numel(qs)
        for s = seeds
            if f == 1   % reused runs, always the real 4000-step ones
                base = sprintf('Qerr%s_snr30_seed%d_it4000_cosine', qs{q}, s);
            else
                base = sprintf('Qerr%s_snr30_seed%d_it%d_cosine_fref25', qs{q}, s, iters);
            end
            T = readtable(fullfile(outDir, ['summary_' base '.csv']), 'TextType', 'char');
            H = readtable(fullfile(outDir, ['history_' base '.csv']), 'TextType', 'char');
            for n = 1:numel(names)
                if n == 1, r = T(startsWith(T.names, 'Tikhonov'), :); best = nan; vm = nan; ve = nan;
                else
                    r = T(strcmp(T.names, names{n}), :);
                    h = H(strcmp(H.config, names{n}), :);
                    best = max(h.snrAll); vm = h.valMC(end); ve = h.valEI(end);
                end
                rows(end+1, :) = {refs{f}, qs{q}, s, names{n}, r.snrAll, r.snrDeep, best, vm, ve}; %#ok<SAGROW>
            end
        end
    end
end
R = cell2table(rows, 'VariableNames', {'fref','Qerr','seed','method','snrAll','snrDeep','bestAll','valMC','valEI'});
writetable(R, fullfile(outDir, 'fref_table.csv'));
disp(R);

get = @(f, q, m, col) arrayfun(@(s) R.(col)(strcmp(R.fref, f) & strcmp(R.Qerr, q) & R.seed == s & strcmp(R.method, m)), seeds);
loss = @(f, q, m) get(f, '+0.00', m, 'snrAll') - get(f, q, m, 'snrAll');
fid = fopen(fullfile(outDir, 'fref_verdicts.txt'), 'w');
pr = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);

h1 = {};
for q = errs
    r = loss('Nyquist', q{1}, 'MC + EI') - loss('25Hz', q{1}, 'MC + EI');
    if all(r > 0) && mean(r) >= 2, v = 'clearly amplified';
    elseif all(r > 0) && mean(r) >= 1, v = 'somewhat amplified';
    else, v = 'no clear amplification'; end
    h1{end+1} = v; %#ok<SAGROW>
    pr('H1 Qerr %s: r for MC+EI = %s (mean %.2f) -> %s\n', q{1}, mat2str(r, 3), mean(r), v);
end
h23 = {};
for q = qs
    d = get('25Hz', q{1}, 'MC + EI', 'snrAll') - get('25Hz', q{1}, 'MC only', 'snrAll'); v2 = rule(d);
    pr('H2 Qerr %s: MC+EI - MC = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), v2);
    d = get('25Hz', q{1}, 'MC + EI', 'snrAll') - get('25Hz', q{1}, 'Tikhonov', 'snrAll'); v3 = rule(d);
    pr('H3 Qerr %s: MC+EI - Tikhonov = %s (mean %.2f) -> %s\n', q{1}, mat2str(d, 3), mean(d), v3);
    if ismember(q{1}, errs), h23 = [h23, {v2, v3}]; end %#ok<AGROW>
end
if all(strcmp(h23, 'advantage retained')), overall = 'A: earlier unfavourable result mainly from the Nyquist reference';
elseif any(strcmp(h1, 'clearly amplified')), overall = 'B: loss amplified, but EI advantage still insufficient';
else, overall = 'C: EI sensitivity to Q is basically not from the reference frequency'; end
pr('Overall: %s\n\nDescriptive, mean over seeds\n', overall);

for m = names
    for f = refs
        pr('%-9s %-8s all  %s | deep %s | loss %s\n', m{1}, f{1}, ...
            sprintf('%7.2f', cellfun(@(q) mean(get(f{1}, q, m{1}, 'snrAll')), qs)), ...
            sprintf('%7.2f', cellfun(@(q) mean(get(f{1}, q, m{1}, 'snrDeep')), qs)), ...
            sprintf('%7.2f', cellfun(@(q) mean(get(f{1}, '+0.00', m{1}, 'snrAll') - get(f{1}, q, m{1}, 'snrAll')), qs)));
    end
    pr('%-9s r (Nyquist loss - 25 Hz loss), Qerr -0.10 / +0.10: %s\n', m{1}, ...
        sprintf('%7.2f', cellfun(@(q) mean(loss('Nyquist', q, m{1}) - loss('25Hz', q, m{1})), errs)));
end
for m = names(2:3)
    pr('%-9s 25Hz best recorded %s | held-out MC %s | held-out EI %s\n', m{1}, ...
        sprintf('%7.2f', cellfun(@(q) mean(get('25Hz', q, m{1}, 'bestAll')), qs)), ...
        sprintf('%8.4f', cellfun(@(q) mean(get('25Hz', q, m{1}, 'valMC')), qs)), ...
        sprintf('%8.4f', cellfun(@(q) mean(get('25Hz', q, m{1}, 'valEI')), qs)));
end
fclose(fid);

function v = rule(d)
if all(d > 0) && mean(d) >= 1, v = 'advantage retained';
elseif all(d > 0), v = 'positive but below 1 dB';
elseif all(d < 0) && mean(d) <= -1, v = 'harmful';
else, v = 'undecided'; end
end
