function run_output_trust_case(outdir,inputpath,stage,kappa)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
input=jsondecode(fileread(inputpath));p=input.public;for f={'vp','h','offset','lambda_grid'},key=f{1};p.(key)=reshape(p.(key),1,[]);end
p.t=p.t(:);assert(~any(isfield(p,{'Q','fpeak','truth','rho','vs'})));assert(~isfolder(outdir));mkdir(outdir);started=tic;
if ~isempty(input.cached_calibration)
 c=input.cached_calibration;est=c.estimate;corrected=c.corrected;cal=c.diagnostics;
 for f={'vp','h','offset','lambda_grid'},key=f{1};corrected.(key)=reshape(corrected.(key),1,[]);end
 corrected.t=corrected.t(:);
else
 policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));[est,corrected,cal]=p1_prepare_candidate(input.calibrationD,p.offset,p,policy);
end
cm=uq_inverse_config(corrected,est);d=prune_dictionary(cm,input.section_offset_m);cov=qcluster_calibration_cov(input.calibrationD,p.offset,corrected,est);Cq=cov.logQ_covariance;dB=cell(1,4);step=.02;
for k=1:4
 cp=cm;cn=cm;cp.Q(k)=cm.Q(k)*exp(step);cn.Q(k)=cm.Q(k)*exp(-step);
 [a,b]=ray_absorption(cm.t,input.section_offset_m,cp);Gp=attenuation_matrix(a+b,cp);[a,b]=ray_absorption(cm.t,input.section_offset_m,cn);Gn=attenuation_matrix(a+b,cn);
 dB{k}=((Gp-Gn)*d.W)./(2*step*d.scale');
end
uq_writejson(fullfile(outdir,'calibration.json'),struct('estimate',est,'corrected',corrected,'diagnostics',cal,'uncertainty',cov));
fprintf('CALIBRATION_READY %s elapsed=%.1f\n',input.case_name,toc(started));
if strcmp(stage,'development'),keys={'QAwareRidge','QAwareCluster','Trust05','Trust10','Trust20'};else,keys={'QAwareRidge','QAwareCluster','OutputTrustNoQ','QAwareOutputTrust'};end
flags=0;maxnormal=0;maxratio=0;total=0;
for id=1:numel(input.datasets)
 data=input.datasets(id);D=data.D;nx=size(D,2);sigma=zeros(1,nx);for j=1:nx,ni=estimate_trace_snr(D(:,j),p.t);sigma(j)=ni.sigma;end
 fprintf('SUPPORT_START %s %s traces=%d\n',input.case_name,data.name,nx);
 initial=reflectivity_refit_prepared(D,repmat({d},1,nx),sigma);pruned=reflectivity_prune_prepared(D,repmat({d},1,nx),sigma,initial);
 flags=flags+initial.optimization_flags;maxnormal=max([maxnormal,initial.normal_residual,pruned.normal_residual]);
 result=struct('Observed',D,'FullCalRefit',initial.gather,'PrunedFullCalRefit',pruned.gather);stats=struct;
 for key=keys
  name=key{1};S=zeros(size(D));di=cell(1,nx);
  for j=1:nx
   if any(strcmp(name,{'QAwareRidge','QAwareCluster'}))
    out=qcluster_refit(D(:,j),d,sigma(j),initial.reflectivity(:,j),initial.support{j},dB,Cq,name);
   else
    kk=kappa;if strcmp(name,'Trust05'),kk=.5;elseif strcmp(name,'Trust10'),kk=1;elseif strcmp(name,'Trust20'),kk=2;end
    useQ=~strcmp(name,'OutputTrustNoQ');out=qtrust_refit(D(:,j),d,sigma(j),initial.reflectivity(:,j),initial.support{j},dB,Cq,kk,useQ);
    maxratio=max(maxratio,out.max_protected_change_sd/max(kk,realmin));
   end
   S(:,j)=out.gather;di{j}=rmfield(out,'gather');maxnormal=max(maxnormal,out.normal_residual);
  end
  result.(name)=S;stats.(name)=di;fprintf('METHOD_DONE %s %s %s rank=%.1f elapsed=%.1f\n',input.case_name,data.name,name,mean(cellfun(@(v)v.rank_after,di)),toc(started));
 end
 result.sigma=sigma;result.diagnostics=stats;result.support=initial.support;result.initial_flags=initial.optimization_flags;
 uq_writejson(fullfile(outdir,[data.name '.json']),result);
 save(fullfile(outdir,[data.name '_prepared.mat']),'d','dB','Cq','D','sigma','initial','pruned','-v7');
 total=total+nx;uq_writejson(fullfile(outdir,'progress.json'),struct('case',input.case_name,'traces',total,'elapsed_seconds',toc(started)));
end
control=qtrust_controls();assert(flags==0&&maxnormal<1e-8&&maxratio<=1+1e-6);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'case',input.case_name,'traces',total,'optimization_flags',flags,'max_normal',maxnormal,'max_guard_ratio',maxratio,'controls',control,'elapsed_seconds',toc(started)));
end
