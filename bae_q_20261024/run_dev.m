function run_dev(which)
% Development run on soft_refit_20261023 scenes (retrospective development data).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
test_bae_math();opts=jsondecode(fileread(fullfile(root,'SPEC.json')));opts.s_grid=reshape(opts.s_grid,1,[]);
soft=fullfile(fileparts(root),'soft_refit_20261023');[cm,ref]=load_dev(soft,which);started=tic;
outdir=fullfile(root,'results','dev');if ~isfolder(outdir),mkdir(outdir);end
for errorpct=reshape(opts.errors_percent,1,[])
 for k=1:numel(ref.datasets)
  data=ref.datasets(k);file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',which,errorpct,data.name));
  if isfile(file),continue;end
  [g,diag]=run_scene(cm,data.D,data.sigma,ref.reference_Q,errorpct,opts);
  names=fieldnames(g);for q=1:numel(names),g.(names{q})=single(g.(names{q}));end
  save(file,'g','diag','-v7');
  fprintf('%s Q%+d %s flags=%d normal=%.2g rank=%d elapsed=%.1fs\n',which,errorpct,data.name,diag.optimization_flags,diag.max_normal_residual,diag.max_error_rank,toc(started));
 end
end
fprintf('DONE %s %.1fs\n',which,toc(started));
end
