%% valmc_signal_posthoc.m -- post hoc, descriptive only; not part of any pre-registered rule.
% Across all 4000-step cosine runs with seeds 1-3: does the held-out MC loss of MC+EI relative to
% MC track the SNR difference? The held-out MC loss uses no labels; the SNR difference does.
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
conds = {'Q correct, stationary', 'Qerr+0.00_snr30_seed%d_it4000_cosine'
         'Qerr +0.2',             'Qerr+0.20_snr30_seed%d_it4000_cosine'
         'Qerr -0.2',             'Qerr-0.20_snr30_seed%d_it4000_cosine'
         'density',               'Qerr+0.00_snr30_seed%d_it4000_cosine_density'
         'amp',                   'Qerr+0.00_snr30_seed%d_it4000_cosine_amp'
         'broad2narrow',          'Qerr+0.00_snr30_seed%d_it4000_cosine_broad2narrow'
         'narrow2broad',          'Qerr+0.00_snr30_seed%d_it4000_cosine_narrow2broad'};
rows = {};
for c = 1:size(conds, 1)
    for s = 1:3
        H = readtable(fullfile(outDir, ['history_' sprintf(conds{c, 2}, s) '.csv']), 'TextType', 'char');
        M = H(strcmp(H.config, 'MC only'), :); E = H(strcmp(H.config, 'MC + EI'), :);
        rows(end+1, :) = {conds{c, 1}, s, E.snrAll(end) - M.snrAll(end), M.valMC(end), E.valMC(end), ...
            E.valMC(end) / M.valMC(end), M.valEI(end), E.valEI(end)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'condition','seed','dSNR_EI_minus_MC','valMC_MC','valMC_EI','valMC_ratio_EI_over_MC','valEI_MC','valEI_EI'});
writetable(T, fullfile(outDir, 'valmc_signal_posthoc.csv'));
disp(T);
fprintf('Corr(dSNR, log ratio) over %d runs: %.3f\n', height(T), corr(T.dSNR_EI_minus_MC, log(T.valMC_ratio_EI_over_MC)));
disp(groupsummary(T, 'condition', 'mean', {'dSNR_EI_minus_MC','valMC_ratio_EI_over_MC'}));
Qok = ~startsWith(T.condition, 'Qerr');
fprintf('Q correct only, %d runs: corr %.3f; ratio > 1 in %d runs, their dSNR %s; ratio <= 1 min dSNR %.2f\n', nnz(Qok), ...
    corr(T.dSNR_EI_minus_MC(Qok), log(T.valMC_ratio_EI_over_MC(Qok))), nnz(Qok & T.valMC_ratio_EI_over_MC > 1), ...
    mat2str(T.dSNR_EI_minus_MC(Qok & T.valMC_ratio_EI_over_MC > 1), 3), min(T.dSNR_EI_minus_MC(Qok & T.valMC_ratio_EI_over_MC <= 1)));
