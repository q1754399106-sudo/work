%% aggregate_multiseed.m -- applies the rules in PLAN_MULTISEED_ZH.md to outputs/*_cosine files.
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
seeds = [1 2 3]; budgets = [4000 12000];
if ~isempty(getenv('EI_AGG_BUDGETS')), budgets = str2num(getenv('EI_AGG_BUDGETS')); end %#ok<ST2NM> % smoke test only
fid = fopen(fullfile(outDir, 'multiseed_verdicts.txt'), 'w');
pr = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
rows = {};
stab = {};
for b = budgets
    for s = seeds
        tag = sprintf('Qerr+0.00_snr30_seed%d_it%d_cosine', s, b);
        S = readtable(fullfile(outDir, ['summary_' tag '.csv']), 'TextType', 'char');
        H = readtable(fullfile(outDir, ['history_' tag '.csv']), 'TextType', 'char');
        cfgs = unique(H.config, 'stable');
        for k = 1:numel(cfgs)
            h = H(strcmp(H.config, cfgs{k}), :);
            last3 = h.snrAll(end-2:end);
            lam = double(~strcmp(cfgs{k}, 'MC only'));           % same weights as the training loss
            [~, iv] = min(h.valMC + lam * h.valEI);
            stab(end+1, :) = {b, s, cfgs{k}, max(last3) - min(last3), min(h.snrAll(ceil(end/2):end)), h.snrAll(end), h.snrDeep(end), h.iter(iv), h.snrAll(iv)}; %#ok<SAGROW>
        end
        tik = S(startsWith(S.names, 'Tikhonov'), :);
        rows(end+1, :) = {b, s, tik.snrAll, tik.snrDeep}; %#ok<SAGROW>
    end
end
T = cell2table(stab, 'VariableNames', {'budget','seed','config','range_last3','min_second_half','snrAll','snrDeep','iter_valbest','snrAll_valbest'});
writetable(T, fullfile(outDir, 'multiseed_table.csv'));
disp(T);
Tik = cell2table(rows, 'VariableNames', {'budget','seed','tikAll','tikDeep'});
disp(Tik);

pr('C1 stability: max range of last 3 records = %.2f dB -> %s\n', max(T.range_last3), yn(max(T.range_last3) <= 1));
get = @(b, cfg, col) arrayfun(@(s) T.(col)(T.budget == b & T.seed == s & strcmp(T.config, cfg)), seeds);
for b = budgets
    d  = get(b, 'MC + EI', 'snrAll')  - get(b, 'MC only', 'snrAll');
    dd = get(b, 'MC + EI', 'snrDeep') - get(b, 'MC only', 'snrDeep');
    if all(d > 0) && mean(d) >= 1, v = 'advantage holds';
    elseif all(d > 0), v = 'positive but below 1 dB';
    else, v = 'undecided (signs differ)'; end
    pr('C2 %d steps: MC+EI - MC = %s (mean %.2f) -> %s | deep %s (mean %.2f)\n', b, mat2str(d, 3), mean(d), v, mat2str(dd, 3), mean(dd));
end
d = get(budgets(1), 'MC + EI + fmax', 'snrAll') - get(budgets(1), 'MC + EI', 'snrAll');
if all(d > 0) && mean(d) >= 0.5, v = 'gain';
elseif all(d < 0) && mean(d) <= -0.5, v = 'harmful';
else, v = 'no gain seen'; end
pr('C3 4000 steps: fmax - MC+EI = %s (mean %.2f) -> %s\n', mat2str(d, 3), mean(d), v);
fclose(fid);

function s = yn(b)
if b, s = 'holds'; else, s = 'fails'; end
end
