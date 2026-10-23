function run_dev_subset(which,errors)
% Same computation as run_dev for a subset of Q errors, so a second MATLAB
% process can fill the grid from the other end. Existing files are skipped.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
test_bae_math();opts=jsondecode(fileread(fullfile(root,'SPEC.json')));opts.s_grid=reshape(opts.s_grid,1,[]);
soft=fullfile(fileparts(root),'soft_refit_20261023');[cm,ref]=load_dev(soft,which);started=tic;outdir=fullfile(root,'results','dev');
for errorpct=reshape(errors,1,[])
 for k=1:numel(ref.datasets)
  data=ref.datasets(k);file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',which,errorpct,data.name));
  if isfile(file),continue;end
  [g,diag]=run_scene(cm,data.D,data.sigma,ref.reference_Q,errorpct,opts);
  names=fieldnames(g);for q=1:numel(names),g.(names{q})=single(g.(names{q}));end
  save(file,'g','diag','-v7');
  fprintf('%s Q%+d %s flags=%d normal=%.2g rank=%d elapsed=%.1fs\n',which,errorpct,data.name,diag.optimization_flags,diag.max_normal_residual,diag.max_error_rank,toc(started));
 end
end
fprintf('DONE subset %s %.1fs\n',which,toc(started));
end
