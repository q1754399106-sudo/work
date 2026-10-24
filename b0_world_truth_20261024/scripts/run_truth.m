function run_truth(u)
% Scoring side, only after the confirmation predictions are sealed: expected Delta MSE of QScan
% versus FullCalRefit at the public Q for model U<u>, 12 noise replicates per SNR level (the
% diagnosed gathers are replicates 1-3), plus the realized errors of the diagnosed gathers.
P=bw_setup();maxNumCompThreads(P.threads);assert(isfile(fullfile(P.root,'predictions','SEALED_SHA256.json')));
in=load(fullfile(P.root,'inputs',sprintf('U%d.mat',u)));tr=load(fullfile(P.root,'truth',sprintf('U%d_truth.mat',u)));
cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);T=tr.target;D0=tr.clean;ET=sum(T(:).^2);snr=[8 11 14 17];rows={};
for level=1:4
 for r=1:12
  seed=1200000000+10000000*u+100000*level+r;D=D0+noise_replica(D0,snr(level),seed);sig=estimate_sigmas(D,cm.t);
  if r<=3,o=in.obs(strcmp({in.obs.name},sprintf('SNR%02d',snr(level)))&[in.obs.replicate]==r);assert(isequal(o.D,D));end
  sel=qscan_select(D,sig,cm,Q0,cm.solver_offset);c=P.grid(sel.best);
  x1=full_refit(D,sig,cm,Q0,cm.solver_offset);xs=full_refit(D,sig,cm,Q0*c,cm.solver_offset);
  rows(end+1,:)={sprintf('U%d',u),sprintf('SNR%02d',snr(level)),r,seed,c,sum((x1-T).^2,'all')/ET,sum((xs-T).^2,'all')/ET}; %#ok<AGROW>
 end
end
writetable(cell2table(rows,'VariableNames',{'model','condition','replicate','seed','qscan_factor','mse_base','mse_qscan'}), ...
 fullfile(P.root,'analysis',sprintf('truth_U%d.csv',u)));
fprintf('DONE truth U%d\n',u);
end
