function [cm,ref,input]=load_dev(softroot,which)
% Development scenes of soft_refit_20261023 (same observations, sigma, config).
input=jsondecode(fileread(fullfile(softroot,'inputs',[which '.json'])));
ref=jsondecode(fileread(fullfile(softroot,'results',which,'reference.json')));
if strcmp(which,'primary')
 cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';
else
 p=input.public;for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end
 p.t=p.t(:);cm=uq_inverse_config(p,input.estimate);
end
end
