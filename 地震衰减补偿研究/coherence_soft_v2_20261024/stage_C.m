function stage_C()
% Registered stage C primary: Q-uncertainty energy in the gain (reference = stage A winner).
ref=cand('A_0.5_16',0);c=ref([]);
for kappa=[.25 .5 1],c(end+1)=cand(sprintf('VQ_k%g',kappa),kappa);end %#ok<AGROW>
run_stage('C',c,ref,@(s,d)q_spread_energy(s,d,[.8 1.25]));
end

function x=cand(label,kappa)
x=struct('label',label,'opts',struct('lambda_min',.5,'lambda_max',16,'tau',0,'kappa',kappa),'cfun',[]);
end
