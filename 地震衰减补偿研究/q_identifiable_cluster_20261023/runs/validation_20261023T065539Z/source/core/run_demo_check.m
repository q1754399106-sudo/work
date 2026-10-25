function run_demo_check(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);
result=demo_optimized(outdir);frozen=jsondecode(fileread(fullfile(root,'inputs','frozen_pruned_profile.json')));
reference=load(fullfile(frozen.development_run,'matlab','cases','OriginalReference_WeakSection_12dB_20261411.mat'));
replay=norm(result.gather-reference.optimized,'fro')/norm(reference.optimized,'fro');assert(replay<1e-8&&result.optimization_flags==0);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'replay',replay,'method',result.selected_method));
end
