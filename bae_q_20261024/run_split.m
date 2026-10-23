function run_split(split,models)
% Frozen run on 'val' (coherence_soft V1-V3) or 'ind' (independent F1-F2).
% Uses only the development-selected s values (analysis/SELECTED.json).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
test_bae_math();sel=jsondecode(fileread(fullfile(root,'analysis','SELECTED.json')));assert(sel.frozen_before_validation);
opts=jsondecode(fileread(fullfile(root,'SPEC.json')));opts.s_grid=reshape(sel.s_values,1,[]);opts.lean=true;
outdir=fullfile(root,'results',split);if ~isfolder(outdir),mkdir(outdir);end
scenes=list_scenes(root,split);started=tic;
for k=1:numel(scenes)
 sc=scenes(k);if nargin>1&&~any(strcmp(sc.model,models)),continue;end
 for errorpct=reshape(opts.errors_percent,1,[])
  file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',sc.model,errorpct,sc.name));if isfile(file),continue;end
  [g,diag]=run_scene(sc.cm,sc.D,sc.sigma,sc.referenceQ,errorpct,opts);
  names=fieldnames(g);for q=1:numel(names),g.(names{q})=single(g.(names{q}));end
  save(file,'g','diag','-v7');
  fprintf('%s Q%+d %s flags=%d normal=%.2g elapsed=%.1fs\n',sc.model,errorpct,sc.name,diag.optimization_flags,diag.max_normal_residual,toc(started));
 end
end
fprintf('DONE %s %.1fs\n',split,toc(started));
end
