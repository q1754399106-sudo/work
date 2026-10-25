function [configs,spec]=p1_models(split_name)
% Return geological model configurations for 'historical', 'development', or 'heldout_new_geology'.
if nargin<1 || isempty(split_name), split_name='historical'; end
root=fileparts(fileparts(mfilename('fullpath')));
spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
mlist=spec.models.(split_name);
configs=cell(1,numel(mlist));
fields={'vp','vs','rho','Q','event_t','fpeak'};
for i=1:numel(mlist)
 c=config_model(); m=mlist(i);
 for k=1:numel(fields), f=fields{k}; c.(f)=reshape(m.(f),1,[]); end
 c.h=diff([0,c.event_t]).*c.vp(1:3)/2; c.model_name=m.name;
 c.nominal_frequency=m.nominal_frequency;
 configs{i}=c;
end
end
