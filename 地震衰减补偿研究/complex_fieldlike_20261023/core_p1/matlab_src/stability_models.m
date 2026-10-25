function [configs,spec]=stability_models()
root=fileparts(fileparts(mfilename('fullpath'))); spec=jsondecode(fileread(fullfile(root,'inputs','stability_spec.json')));
configs=cell(1,numel(spec.models)+1); configs{1}=config_model(); configs{1}.model_name='OriginalControl';
fields={'vp','vs','rho','Q','event_t','fpeak'};
for i=1:numel(spec.models)
 c=config_model(); m=spec.models(i);
 for k=1:numel(fields), f=fields{k}; c.(f)=reshape(m.(f),1,[]); end
 c.h=diff([0,c.event_t]).*c.vp(1:3)/2; c.model_name=m.name; configs{i+1}=c;
end
end
