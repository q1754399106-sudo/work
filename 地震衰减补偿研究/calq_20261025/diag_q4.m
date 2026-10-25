function diag_q4()
% Why is the deepest calibrated Q off? Identifiability (objective profile in Q4),
% noise-free bias and noise-replicate spread of the frozen layered calibration.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));o=jsondecode(fileread(fullfile(root,'inputs','uq_spec.json')));
out=struct();started=tic;
for v={'V3','V5'}
 in=load(fullfile(root,'dev',[v{1} '.mat']));TT=load(fullfile(root,'dev',[v{1} '_truth.mat']),'g');model=gen_complex_truth(TT.g);
 for k=1:2,dv=norm(model.datasets(k).D-in.obs(k).D,'fro')/norm(in.obs(k).D,'fro');assert(dv<1e-10,'regeneration differs');end
 p=model.public;qtrue=model.truth.calibration_Q;
 [est,corrected]=p1_prepare_candidate(model.calibrationD,p.offset,p,policy);
 fprintf('%s archived Q %s | recomputed %s | true coarse %s\n',v{1},mat2str(in.config.Q,4),mat2str(est.layer_Q,4),mat2str(qtrue,4));
 % Objective profile in Q4 (other parameters at the estimate), with and without the smoothness prior.
 q4=exp(linspace(log(15),log(200),40));Lp=zeros(size(q4));Ld=Lp;
 for i=1:numel(q4),qs=est.layer_Q;qs(4)=q4(i);[Lp(i),Ld(i)]=layered_loss(model.calibrationD,corrected,qs,est.estimated_frequency,o);end
 % Noise-free calibration and noise replicates (new seeds).
 [e0,~]=p1_prepare_candidate(model.calibration_clean,p.offset,p,policy);
 reps=zeros(8,4);clean=model.calibration_clean;
 for r=1:8
  rng(2026102500+r,'twister');D=clean;
  for j=1:11,n=randn(p.nt,1);D(:,j)=clean(:,j)+n/sqrt(mean(n.^2))*sqrt(mean(clean(:,j).^2))*10^(-18/20);end
  e=p1_prepare_candidate(D,p.offset,p,policy);reps(r,:)=e.layer_Q;
 end
 lr=log(reps./qtrue);
 fprintf('%s noise-free estimate %s\n',v{1},mat2str(e0.layer_Q,4));
 fprintf('%s replicates ln(Q/Qtrue): mean %s  std %s\n',v{1},mat2str(mean(lr,1),3),mat2str(std(lr,0,1),3));
 [~,ip]=min(Lp);[~,id]=min(Ld);
 fprintf('%s Q4 profile: with prior min at %.1f (range of loss %.3g), data-only min at %.1f (range %.3g) | %.0fs\n',v{1},q4(ip),max(Lp)-min(Lp),q4(id),max(Ld)-min(Ld),toc(started));
 out.(v{1})=struct('archived_Q',in.config.Q,'recomputed_Q',est.layer_Q,'true_coarse_Q',qtrue,'noise_free_Q',e0.layer_Q, ...
  'replicate_Q',reps,'replicate_log_bias',mean(lr,1),'replicate_log_std',std(lr,0,1),'q4_grid',q4,'loss_with_prior',Lp,'loss_data_only',Ld, ...
  'calibration_event_times',cumsum(2*p.h./p.vp(1:3)));
end
uq_writejson(fullfile(root,'analysis','diag_q4.json'),out);
end

function [withprior,dataonly]=layered_loss(D,public,qs,fp,o)
% Same loss as p1_layered_calibrate (training traces), split into prior and data parts.
train=reshape(o.calibration_trace_indices,1,[]);times=cumsum(2*public.h./public.vp(1:3));mask=false(public.nt,1);
for t=times,mask=mask|abs(public.t-t)<=o.event_half_window_s;end
sigma=zeros(1,size(D,2));for j=1:size(D,2),ni=estimate_trace_snr(D(:,j),public.t);sigma(j)=ni.sigma;end
dataonly=0;
for j=train
 B=p1_layered_calibration_templates(public,public.offset(j),qs,fp);X=B(mask,:)/sigma(j);y=D(mask,j)/sigma(j);w=ones(size(y));
 for it=1:o.irls_iterations,G=X'*(w.*X);coef=(G+1e-8*max(trace(G)/6,realmin)*eye(6))\(X'*(w.*y));r=y-X*coef;w=min(1,o.huber_delta./max(abs(r),realmin));end
 r=abs(y-X*coef);s=r<=o.huber_delta;dataonly=dataonly+sum(.5*r(s).^2)+sum(o.huber_delta*(r(~s)-.5*o.huber_delta));
end
prior=o.source_prior_weight*(log(fp/public.nominal_frequency)/o.source_prior_log_scale)^2+.5*sum(diff(log(qs)).^2);
withprior=dataonly+prior;
end
