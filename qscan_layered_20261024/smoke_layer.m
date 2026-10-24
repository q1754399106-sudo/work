function smoke_layer()
% Smoke test on 22 development traces (F1 Moderate, Q error 0): runs, timing, nesting. Not a result.
root=fileparts(mfilename('fullpath'));up=fullfile(root,'upstream');addpath(fullfile(root,'core'),fullfile(up,'bae_core'),fullfile(up,'bae_source'));maxNumCompThreads(4);
spec=jsondecode(fileread(fullfile(root,'SPEC.json')));in=load(fullfile(up,'independent','F1.mat'));cm=in.config;cm.t=cm.t(:);
tic;r=layer_scan(in.obs(1).D(:,1:22),in.obs(1).sigma(1:22),cm,0,spec);fprintf('time %.1fs evals %d c0 %.3f\n',toc,r.evaluations,r.c0);
for v=r.variants,fprintf('%s %s\n',v.name,num2str(v.factors,'%.3f '));end
end
