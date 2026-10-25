function run_truth_oct(u)
% scripts/run_truth.m on Octave, after the seal. Same seeds, replicates, estimators and columns.
% Replicates 1-3 are the stored diagnosed gathers: the regenerated gather is compared with the
% stored one by tolerance (randn shim and BLAS can differ from MATLAB in the last bits) and the
% stored gather is used, so the truth refers to exactly the data the predictions saw.
P=oct_setup();assert(isfile(fullfile(P.root,'predictions','SEALED_SHA256.json')));
in=load(fullfile(P.root,'inputs',sprintf('U%d.mat',u)));tr=load(fullfile(P.root,'truth',sprintf('U%d_truth.mat',u)));
cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);T=tr.target;D0=tr.clean;ET=sum(T(:).^2);snr=[8 11 14 17];rows=[];
for level=1:4
 for r=1:12
  seed=1200000000+10000000*u+100000*level+r;D=D0+noise_replica(D0,snr(level),seed);dev=NaN;
  if r<=3
   o=in.obs(strcmp({in.obs.name},sprintf('SNR%02d',snr(level)))&[in.obs.replicate]==r);
   dev=max(abs(o.D(:)-D(:)))/max(abs(o.D(:)));assert(dev<=1e-12,'regenerated gather differs: %g',dev);D=o.D;
  end
  sig=estimate_sigmas(D,cm.t);sel=qscan_select(D,sig,cm,Q0,cm.solver_offset);c=P.grid(sel.best);
  x1=full_refit(D,sig,cm,Q0,cm.solver_offset);xs=full_refit(D,sig,cm,Q0*c,cm.solver_offset);
  row=struct('model',sprintf('U%d',u),'condition',sprintf('SNR%02d',snr(level)),'replicate',r,'seed',seed,'qscan_factor',c, ...
   'mse_base',sum((x1(:)-T(:)).^2)/ET,'mse_qscan',sum((xs(:)-T(:)).^2)/ET,'regen_rel_dev',dev);rows=[rows;row]; %#ok<AGROW>
  fprintf('U%d SNR%02d r%02d c=%.4f base %.4f qscan %.4f dev %.1e\n',u,snr(level),r,c,row.mse_base,row.mse_qscan,dev);
 end
end
write_rows_csv(fullfile(P.root,'analysis',sprintf('truth_U%d.csv',u)),rows);
fprintf('DONE truth U%d\n',u);
end
