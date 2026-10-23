function run_qcluster_case(outdir,inputpath)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
input=jsondecode(fileread(inputpath));p=input.public;
for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end
p.t=p.t(:);assert(~any(isfield(p,{'Q','fpeak','truth','rho','vs'})));
assert(~isfolder(outdir));mkdir(outdir);started=tic;
if ~isempty(input.cached_calibration)
 cal0=input.cached_calibration;est=cal0.estimate;corrected=cal0.corrected;
 for field={'vp','h','offset','lambda_grid'},f=field{1};corrected.(f)=reshape(corrected.(f),1,[]);end
 corrected.t=corrected.t(:);cal=cal0.diagnostics;
else
 policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
 [est,corrected,cal]=p1_prepare_candidate(input.calibrationD,p.offset,p,policy);
end
cov=qcluster_calibration_cov(input.calibrationD,p.offset,corrected,est);Cq=cov.logQ_covariance;
if isfield(input,'q_point_scale')
 est.layer_Q=reshape(est.layer_Q,1,[])*input.q_point_scale;est.effective_Q=est.effective_Q*input.q_point_scale;
end
cm=uq_inverse_config(corrected,est);d=prune_dictionary(cm,input.section_offset_m);
dB=cell(1,4);step=.02;
for k=1:4
 cp=cm;cn=cm;cp.Q(k)=cm.Q(k)*exp(step);cn.Q(k)=cm.Q(k)*exp(-step);
 [a,b]=ray_absorption(cm.t,input.section_offset_m,cp);Gp=attenuation_matrix(a+b,cp);
 [a,b]=ray_absorption(cm.t,input.section_offset_m,cn);Gn=attenuation_matrix(a+b,cn);
 dB{k}=((Gp-Gn)*d.W)./(2*step*d.scale');
end
uq_writejson(fullfile(outdir,'calibration.json'),struct('estimate',est,'corrected',corrected,'diagnostics',cal,'uncertainty',cov));
fprintf('CALIBRATION_READY %s elapsed=%.1f Q=[%s]\n',input.case_name,toc(started),num2str(cm.Q));
methods={'GlobalTSVD','ClusterModesNoQ','QAwareRidge','QAwareCluster'};flags=0;maxnormal=0;total=0;
for id=1:numel(input.datasets)
 data=input.datasets(id);D=data.D;nx=size(D,2);sigma=zeros(1,nx);
 for j=1:nx,ni=estimate_trace_snr(D(:,j),p.t);sigma(j)=ni.sigma;end
 fprintf('SUPPORT_START %s %s traces=%d\n',input.case_name,data.name,nx);
 initial=reflectivity_refit_prepared(D,repmat({d},1,nx),sigma);
 pruned=reflectivity_prune_prepared(D,repmat({d},1,nx),sigma,initial);
 flags=flags+initial.optimization_flags;maxnormal=max([maxnormal,initial.normal_residual,pruned.normal_residual]);
 result=struct('Observed',D,'FullCalRefit',initial.gather,'PrunedFullCalRefit',pruned.gather);
 diagnostics=struct;
 for method=methods
  key=method{1};S=zeros(size(D));stats=cell(1,nx);
  for j=1:nx
   out=qcluster_refit(D(:,j),d,sigma(j),initial.reflectivity(:,j),initial.support{j},dB,Cq,key);
   S(:,j)=out.gather;stats{j}=rmfield(out,'gather');maxnormal=max(maxnormal,out.normal_residual);
  end
  result.(key)=S;diagnostics.(key)=stats;
  fprintf('METHOD_DONE %s %s %s rank_mean=%.2f elapsed=%.1f\n',input.case_name,data.name,key,mean(cellfun(@(x)x.rank_after,stats)),toc(started));
 end
 result.sigma=sigma;result.diagnostics=diagnostics;result.support=initial.support;result.initial_flags=initial.optimization_flags;result.initial_normal=initial.normal_residual;
 uq_writejson(fullfile(outdir,[data.name '.json']),result);total=total+nx;
 uq_writejson(fullfile(outdir,'progress.json'),struct('case',input.case_name,'dataset',data.name,'traces',total,'elapsed_seconds',toc(started)));
end
% Limiting-case check on real prepared matrices: zero Q covariance must exactly remove whitening.
j=find(vecnorm(D)>0,1);a=qcluster_refit(D(:,j),d,sigma(j),initial.reflectivity(:,j),initial.support{j},dB,zeros(4),'QAwareCluster');
b=qcluster_refit(D(:,j),d,sigma(j),initial.reflectivity(:,j),initial.support{j},dB,zeros(4),'ClusterModesNoQ');
zero_limit=norm(a.gather-b.gather)/max(norm(b.gather),realmin);assert(zero_limit<1e-10);
c=qcluster_refit(D(:,j),d,sigma(j),initial.reflectivity(:,j),initial.support{j},dB,zeros(4),'QAwareRidge');
refit_limit=norm(c.gather-initial.gather(:,j))/max(norm(initial.gather(:,j)),realmin);assert(refit_limit<1e-6);
z=qcluster_refit(zeros(p.nt,1),d,0,zeros(p.nt,1),[],dB,Cq,'QAwareCluster');assert(norm(z.gather)==0);
controls=qcluster_controls();
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'case',input.case_name,'traces',total,'optimization_flags',flags,'max_normal',maxnormal,'zero_Q_covariance_limit',zero_limit,'zero_Q_ridge_refit_limit',refit_limit,'zero_input_exact',norm(z.gather)==0,'synthetic_controls',controls,'elapsed_seconds',toc(started)));
end
