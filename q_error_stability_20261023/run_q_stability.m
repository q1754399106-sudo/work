function run_q_stability(which)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));started=tic;
if strcmp(which,'primary')
 input=jsondecode(fileread(fullfile(root,'inputs','primary.json')));c=input.config;c.t=c.t(:);c.vp=c.vp(:)';c.h=c.h(:)';c.Q=c.Q(:)';
 R=input.reflectivity;t=c.t;W=ricker(c);target=W*R;[aa,bb]=ray_absorption(t,450,c);G=attenuation_matrix(aa+bb,c);clean=G*target;
 datasets=struct('name',{},'D',{},'clean',{},'noise',{},'target',{},'sigma',{});
 for k=1:numel(input.noise)
  item=input.noise(k);N=item.unit_noise;N=N*(norm(clean,'fro')/(10^(item.snr_db/20)*norm(N,'fro')));D=clean+N;
  sigma=zeros(1,size(D,2));for j=1:size(D,2),ni=estimate_trace_snr(D(:,j),t);sigma(j)=ni.sigma;end
  datasets(k)=struct('name',item.name,'D',D,'clean',clean,'noise',N,'target',target,'sigma',sigma);
 end
 cm=c;reference_Q=c.Q;scenarios=input.scenarios;
elseif strcmp(which,'complex')
 input=jsondecode(fileread(fullfile(root,'inputs','complex.json')));p=input.public;
 for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end;p.t=p.t(:);
 cm=uq_inverse_config(p,input.estimate);reference_Q=cm.Q;scenarios=input.scenarios;t=cm.t;
 datasets=input.datasets;
 for k=1:numel(datasets)
  D=datasets(k).D;sigma=zeros(1,size(D,2));for j=1:size(D,2),ni=estimate_trace_snr(D(:,j),t);sigma(j)=ni.sigma;end
  datasets(k).sigma=sigma;
 end
else,error('Unknown protocol');end
outdir=fullfile(root,['matlab_' which]);assert(~isfolder(outdir));mkdir(outdir);
uq_writejson(fullfile(outdir,'reference.json'),struct('t',t,'datasets',datasets,'reference_Q',reference_Q,'protocol',which));
flags=0;maxnormal=0;done=0;
for is=1:numel(scenarios)
 scenario=scenarios(is);used=cm;used.Q=reference_Q.*reshape(scenario.factors,1,[]);
 op=fixed_v1_operators(used,repmat(450,1,11));dict=prune_dictionary(used,450);dicts=repmat({dict},1,11);
 records=cell(1,numel(datasets));
 for k=1:numel(datasets)
  data=datasets(k);D=data.D;sigma=data.sigma;
  [fixed,unpruned,optimized,diag]=apply_all(D,op,dicts,sigma);
  if strcmp(which,'primary')
   [fc,uc,oc,dc]=apply_all(data.clean,op,dicts,sigma);
  else,fc=[];uc=[];oc=[];dc=struct('flags',0,'normal',0);end
  records{k}=struct('name',data.name,'fixed',fixed,'unpruned',unpruned,'optimized',optimized,'fixed_clean',fc,'unpruned_clean',uc,'optimized_clean',oc,'diagnostics',diag,'clean_diagnostics',dc);
  flags=flags+diag.flags+dc.flags;maxnormal=max([maxnormal,diag.normal,dc.normal]);done=done+1;
 end
uq_writejson(fullfile(outdir,[scenario.name '.json']),struct('scenario',scenario,'used_Q',used.Q,'reference_Q',reference_Q,'records',{records}));
uq_writejson(fullfile(outdir,'progress.json'),struct('protocol',which,'scenarios',is,'expected',numel(scenarios),'dataset_cases',done,'flags',flags,'elapsed_seconds',toc(started)));
fprintf('%s scenario %d/%d %s cases=%d flags=%d elapsed=%.1fs\n',which,is,numel(scenarios),scenario.name,done,flags,toc(started));
end
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'protocol',which,'scenarios',numel(scenarios),'dataset_cases',done,'optimization_flags',flags,'maximum_normal_residual',maxnormal,'elapsed_seconds',toc(started)));
end
function W=ricker(c)
t=c.t(:);z=pi*c.fpeak*(t-t');W=(1-2*z.^2).*exp(-z.^2);
end
function [fixed,unpruned,optimized,diag]=apply_all(D,op,dicts,sigma)
[fixed,norm0,~]=fixed_v1_apply(D,op);initial=reflectivity_refit_prepared(D,dicts,sigma);out=reflectivity_prune_prepared(D,dicts,sigma,initial);
unpruned=initial.gather;optimized=out.gather;
diag=struct('flags',out.optimization_flags,'normal',max(norm0,out.normal_residual),'support_before',out.support_before_pruning,'support_after',out.support_after,'removed_coefficients',out.removed_coefficients);
end
