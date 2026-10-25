function run_confirm_oct(worker,workers)
% scripts/run_confirm.m on Octave: same gathers, ids, seeds, diagnose2 call and saved outputs;
% only the table output is replaced by write_rows_csv (Octave has no struct2table/writetable).
P=oct_setup();assert(isfile(fullfile(P.root,'logs','PRE_CONFIRM_HASHES.json')));
assert(isfile(fullfile(P.root,'logs','RESUME_ADDENDUM_SHA256.json')),'Freeze RESUME_ZH.md first.');
outs=fullfile(P.root,'cache','outputs');if ~isfolder(outs),mkdir(outs);end
pred=fullfile(P.root,'predictions');if ~isfolder(pred),mkdir(pred);end
csv=fullfile(pred,sprintf('confirm_w%d.csv',worker));assert(~isfile(csv),'Worker output exists; refusing to overwrite.');
rows=[];started=tic;id=0;
for u=1:4
 in=load(fullfile(P.root,'inputs',sprintf('U%d.mat',u)));assert(~isfield(in,'target')&&~isfield(in,'clean'));
 cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);
 for k=1:numel(in.obs)
  id=id+1;if mod(id-1,workers)~=worker-1,continue;end
  o=in.obs(k);[r,x1,xs]=diagnose2(P,o.D,cm,Q0,cm.solver_offset,1400000000+1000*id);x1=single(x1);xs=single(xs);
  save('-v7',fullfile(outs,sprintf('c%03d.mat',id)),'x1','xs');
  r.id=id;r.model=sprintf('U%d',u);r.condition=o.name;r.replicate=o.replicate;r.seed=o.seed;rows=[rows;r]; %#ok<AGROW>
  write_rows_csv(csv,rows);
  fprintf('c%03d U%d %s r%d: c=%.4f B0 %+.3f B0M %+.3f B0C %+.3f gain %.4f (%.0fs, total %.0fs)\n',id,u,o.name,o.replicate, ...
   r.qscan_factor,r.B0_rel_delta_mse,r.B0M_rel_delta_mse,r.B0C_rel_delta_mse,r.gate_gain,r.seconds,toc(started));
 end
end
fprintf('DONE confirm worker %d %.0fs\n',worker,toc(started));
end
