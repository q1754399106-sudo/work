function run_qscan(split)
% Experiment 1b (registered in PLAN_ZH.md before running): QScanRefit on any split.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
opts=jsondecode(fileread(fullfile(root,'SPEC.json')));outdir=fullfile(root,'results',[split '_qscan']);
if ~isfolder(outdir),mkdir(outdir);end;scenes=list_scenes(root,split);started=tic;
for k=1:numel(scenes)
 sc=scenes(k);offset=450;if isfield(sc.cm,'solver_offset'),offset=sc.cm.solver_offset;end
 for errorpct=reshape(opts.errors_percent,1,[])
  file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',sc.model,errorpct,sc.name));if isfile(file),continue;end
  used=sc.cm;used.Q=reshape(sc.referenceQ,1,[])*(1+errorpct/100);
  [S,info]=qscan_refit(sc.D,sc.sigma,used,offset);g=struct('QScanRefit',single(S));save(file,'g','info','-v7');
  fprintf('%s Q%+d %s factor=%.3f (ideal %.3f) elapsed=%.1fs\n',sc.model,errorpct,sc.name,info.selected_factor,1/(1+errorpct/100),toc(started));
 end
end
fprintf('DONE qscan %s %.1fs\n',split,toc(started));
end
