function out=demo_soft()
% Reproducible whole-gather demonstration; target is not passed to inference.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));
input=jsondecode(fileread(fullfile(root,'inputs','primary.json')));
reference=jsondecode(fileread(fullfile(root,'results','primary','reference.json')));
cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';
data=reference.datasets(2);
out=soft_refit_compensate(data.D,cm,450,reshape(data.sigma,1,[]));
save(fullfile(root,'demo_result.mat'),'out','-v7');
fprintf('Demo complete: %d samples x %d traces.\n',size(out.gather,1),size(out.gather,2));
end
