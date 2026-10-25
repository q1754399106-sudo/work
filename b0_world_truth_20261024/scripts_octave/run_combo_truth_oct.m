function run_combo_truth_oct(u)
% Descriptive only (RESUME_ZH.md), after the seal: combination-level truth on U<u>. Same 12
% replicates as run_truth_oct; "Standard" coherence output at Q0 (B) and at the QScan factor
% (D = QScan + Standard), via scripts/coherence_reference.m (round-3 Standard, lambda 0.5..16).
P=oct_setup();assert(isfile(fullfile(P.root,'predictions','SEALED_SHA256.json')));
in=load(fullfile(P.root,'inputs',sprintf('U%d.mat',u)));tr=load(fullfile(P.root,'truth',sprintf('U%d_truth.mat',u)));
cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);T=tr.target;D0=tr.clean;ET=sum(T(:).^2);snr=[8 11 14 17];rows=[];off=cm.solver_offset;
for level=1:4
 for r=1:12
  seed=1200000000+10000000*u+100000*level+r;D=D0+noise_replica(D0,snr(level),seed);
  if r<=3,o=in.obs(strcmp({in.obs.name},sprintf('SNR%02d',snr(level)))&[in.obs.replicate]==r);D=o.D;end
  sig=estimate_sigmas(D,cm.t);sel=qscan_select(D,sig,cm,Q0,off);c=P.grid(sel.best);
  [XB,x1]=coherence_reference(P,D,sig,cm,Q0,off);[XD,xs]=coherence_reference(P,D,sig,cm,Q0*c,off);
  e=@(X)sum((X(:)-T(:)).^2)/ET;
  row=struct('model',sprintf('U%d',u),'condition',sprintf('SNR%02d',snr(level)),'replicate',r,'seed',seed,'qscan_factor',c, ...
   'mse_base',e(x1),'mse_qscan',e(xs),'mse_standard',e(XB),'mse_qscan_standard',e(XD));rows=[rows;row]; %#ok<AGROW>
  fprintf('U%d SNR%02d r%02d c=%.4f FCR %.4f/%.4f STD %.4f/%.4f\n',u,snr(level),r,c,row.mse_base,row.mse_qscan,row.mse_standard,row.mse_qscan_standard);
 end
end
write_rows_csv(fullfile(P.root,'analysis',sprintf('combo_truth_U%d.csv',u)),rows);
fprintf('DONE combo truth U%d\n',u);
end
