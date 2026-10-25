function run_confirm(worker,workers)
% Confirmation worker: the 48 diagnosed gathers of U1-U4 (4 SNR levels x replicates 1-3), read
% from the solver-side inputs only. Writes truth-free predictions and the two outputs per gather.
P=bw_setup();maxNumCompThreads(P.threads);assert(isfile(fullfile(P.root,'logs','PRE_CONFIRM_HASHES.json')));
outs=fullfile(P.root,'cache','outputs');if ~isfolder(outs),mkdir(outs);end
csv=fullfile(P.root,'predictions',sprintf('confirm_w%d.csv',worker));assert(~isfile(csv),'Worker output exists; refusing to overwrite.');
rows=[];started=tic;id=0;
for u=1:4
 in=load(fullfile(P.root,'inputs',sprintf('U%d.mat',u)));assert(~isfield(in,'target')&&~isfield(in,'clean'));
 cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);
 for k=1:numel(in.obs)
  id=id+1;if mod(id-1,workers)~=worker-1,continue;end
  o=in.obs(k);[r,x1,xs]=diagnose2(P,o.D,cm,Q0,cm.solver_offset,1400000000+1000*id);x1=single(x1);xs=single(xs); %#ok<NASGU>
  save(fullfile(outs,sprintf('c%03d.mat',id)),'x1','xs','-v7');
  r.id=id;r.model=sprintf('U%d',u);r.condition=o.name;r.replicate=o.replicate;r.seed=o.seed;rows=[rows;struct2table(r)]; %#ok<AGROW>
  writetable(rows,csv);
  fprintf('c%03d U%d %s r%d: c=%.4f B0 %+.3f B0M %+.3f B0C %+.3f gain %.4f (%.0fs, total %.0fs)\n',id,u,o.name,o.replicate, ...
   r.qscan_factor,r.B0_rel_delta_mse,r.B0M_rel_delta_mse,r.B0C_rel_delta_mse,r.gate_gain,r.seconds,toc(started));
 end
end
fprintf('DONE confirm worker %d %.0fs\n',worker,toc(started));
end
