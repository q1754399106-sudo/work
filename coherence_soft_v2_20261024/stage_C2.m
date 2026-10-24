function stage_C2()
% Post-interim-check candidates (added after stage C showed losses at correct/high Q):
% coherence-weighted VQ (two-sided spread, cached from stage C) ...
ref=cand('A_0.5_16',0,false);c=ref([]);
for kappa=[.25 .5 1 2],c(end+1)=cand(sprintf('VQcoh_k%g',kappa),kappa,true);end %#ok<AGROW>
c(end+1)=cand('VQ_k0.25',.25,false);
run_stage('C',c,ref,@(s,d)q_spread_energy(s,d,[.8 1.25]));
movefile_safe('C','C2a');
% ... and one-sided spread (Q x 1.25 only: risk of over-compensation).
c=ref([]);for kappa=[.25 .5],c(end+1)=cand(sprintf('VQ1s_k%g',kappa),kappa,false);c(end+1)=cand(sprintf('VQ1scoh_k%g',kappa),kappa,true);end %#ok<AGROW>
run_stage('C1s',c,ref,@(s,d)q_spread_energy(s,d,1.25));
end

function x=cand(label,kappa,coh)
x=struct('label',label,'opts',struct('lambda_min',.5,'lambda_max',16,'tau',0,'kappa',kappa,'vq_coherence',coh),'cfun',[]);
end

function movefile_safe(a,b)
root=fileparts(mfilename('fullpath'));
for f={'metrics.csv','summary.json'}
 movefile(fullfile(root,'analysis',sprintf('stage_%s_%s',a,f{1})),fullfile(root,'analysis',sprintf('stage_%s_%s',b,f{1})),'f');
end
end
