function run_dense()
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));maxNumCompThreads(4);started=tic;flags=0;maxnormal=0;total=0;
for protocol={'primary','complex'}
 which=protocol{1};input=jsondecode(fileread(fullfile(root,'inputs',[which '.json'])));
 if strcmp(which,'primary')
  cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';t=cm.t;z=pi*cm.fpeak*(t-t');W=(1-2*z.^2).*exp(-z.^2);target=W*input.reflectivity;
  [a,b]=ray_absorption(t,450,cm);clean=attenuation_matrix(a+b,cm)*target;
  datasets=struct('name',{},'D',{},'target',{},'sigma',{});
  for k=1:numel(input.noise)
   ns=input.noise(k);N=ns.unit_noise;N=N*(norm(clean,'fro')/(10^(ns.snr_db/20)*norm(N,'fro')));D=clean+N;sigma=get_sigma(D,t);
   datasets(k)=struct('name',ns.name,'D',D,'target',target,'sigma',sigma);
  end
 else
  p=input.public;for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end;p.t=p.t(:);cm=uq_inverse_config(p,input.estimate);t=cm.t;datasets=input.datasets;
  for k=1:numel(datasets),datasets(k).sigma=get_sigma(datasets(k).D,t);end
 end
 referenceQ=cm.Q;outdir=fullfile(root,'results',which);assert(~isfolder(outdir));mkdir(outdir);
 uq_writejson(fullfile(outdir,'reference.json'),struct('t',t,'x',input.x,'datasets',datasets,'reference_Q',referenceQ,'protocol',which));
 for ie=1:numel(input.errors_percent)
  errorpct=input.errors_percent(ie);used=cm;used.Q=referenceQ*(1+errorpct/100);op=fixed_v1_operators(used,repmat(450,1,11));dict=prune_dictionary(used,450);dicts=repmat({dict},1,11);
  records=cell(1,numel(datasets));
  for k=1:numel(datasets)
   data=datasets(k);D=data.D;assert(size(D,2)==99);fixed=zeros(size(D));unpruned=fixed;optimized=fixed;diagnostics=cell(1,9);
   for block=1:9
    ids=(block-1)*11+(1:11);[fixed(:,ids),n0,~]=fixed_v1_apply(D(:,ids),op);initial=reflectivity_refit_prepared(D(:,ids),dicts,data.sigma(ids));out=reflectivity_prune_prepared(D(:,ids),dicts,data.sigma(ids),initial);
    unpruned(:,ids)=initial.gather;optimized(:,ids)=out.gather;flags=flags+out.optimization_flags;maxnormal=max([maxnormal,n0,out.normal_residual]);total=total+1;
    diagnostics{block}=struct('flags',out.optimization_flags,'normal',out.normal_residual,'support_before',out.support_before_pruning,'support_after',out.support_after);
    uq_writejson(fullfile(root,'progress.json'),struct('protocol',which,'error_percent',errorpct,'condition',data.name,'block',block,'blocks_done',total,'blocks_expected',180,'flags',flags,'elapsed_seconds',toc(started)));
    fprintf('%s Q%+d %s block %d/9 total=%d/180 flags=%d elapsed=%.1fs\n',which,errorpct,data.name,block,total,flags,toc(started));
   end
   records{k}=struct('name',data.name,'fixed',fixed,'unpruned',unpruned,'optimized',optimized,'diagnostics',{diagnostics});
  end
  uq_writejson(fullfile(outdir,sprintf('Q_%+d.json',errorpct)),struct('error_percent',errorpct,'used_Q',used.Q,'records',{records}));
 end
 uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'traces',99,'errors_percent',input.errors_percent));
end
uq_writejson(fullfile(root,'checks.json'),struct('complete',true,'blocks',total,'expected',180,'optimization_flags',flags,'maximum_normal_residual',maxnormal,'elapsed_seconds',toc(started)));
end
function sigma=get_sigma(D,t)
sigma=zeros(1,size(D,2));for j=1:size(D,2),if all(D(:,j)==0),sigma(j)=realmin;else,ni=estimate_trace_snr(D(:,j),t);sigma(j)=ni.sigma;end;end
end
