function verify_public_entry()
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));maxNumCompThreads(4);errors=zeros(1,2);k=0;
for protocol={'primary','complex'}
 which=protocol{1};input=jsondecode(fileread(fullfile(root,'inputs',[which '.json'])));
 ref=jsondecode(fileread(fullfile(root,'results',which,'reference.json')));
 expected=jsondecode(fileread(fullfile(root,'results',which,'Q_+0.json')));
 if strcmp(which,'primary')
  cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';
 else
  p=input.public;for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end
  p.t=p.t(:);cm=uq_inverse_config(p,input.estimate);
 end
 cm.Q=reshape(ref.reference_Q,1,[]);data=ref.datasets(2);
 out=soft_refit_compensate(data.D,cm,450,reshape(data.sigma,1,[]));k=k+1;
 errors(k)=norm(out.gather-expected.records(2).soft,'fro')/norm(expected.records(2).soft,'fro');
 assert(errors(k)<1e-12&&out.optimization_flags==0);
 fprintf('Whole-gather API %s relative difference %.3g\n',which,errors(k));
end
uq_writejson(fullfile(root,'API_CHECK.json'),struct('complete',true,'models',2,'traces_each',99,'maximum_relative_difference',max(errors)));
end
