function run_complete_job(jobfile,outdir)
% Experiment harness only: use the submitted seven numerical branches unchanged.
root=fileparts(mfilename('fullpath')); oldpath=path; oldrng=rng;
cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir),'Fresh output directory required'); mkdir(outdir); start=tic;
job=jsondecode(fileread(jobfile)); uq_writejson(fullfile(outdir,'executed_job.json'),job);
c=config_model(); m=job.model;
for key={'vp','vs','rho','Q','event_t','fpeak'}
 k=key{1}; c.(k)=reshape(m.(k),1,[]);
end
c.h=diff([0,c.event_t]).*c.vp(1:3)/2; c.model_name=m.name; c.nominal_frequency=m.nominal_frequency;
public=uq_public_config(c,c.nominal_frequency); trueModels=prepare_models_v2(c.offset,c);
[~,~,~,calScene]=uq_case('SmoothAVO',c,trueModels);
count=0; obs_count=0; maxnormal=0; flags=0; rejected_inputs=0;
previous=[]; audit=struct([]);
for ic=1:numel(job.calibration_seeds)
 cseed=job.calibration_seeds(ic);
 [calD,~]=uq_observation(calScene,c,trueModels,1:11,cseed,18);
 est=p1_calibrate_ensemble(calD,c.offset,public); cm=uq_inverse_config(public,est);
 invmodels=prepare_models_v2(c.offset,cm); R=avo_projector_v2(invmodels);
 delta=0; if ~isempty(previous), delta=norm(invmodels{6}.G-previous{6}.G,'fro'); end
 previous=invmodels;
 uq_writejson(fullfile(outdir,sprintf('calibration_%d.json',cseed)),est);
 save(fullfile(outdir,sprintf('calibration_input_%d.mat',cseed)),'calD','public','-v7');
 record=struct('calibration_seed',cseed,'operator_change',delta,'Q',est.effective_Q, ...
  'frequency',est.estimated_frequency,'q_boundary',est.q_boundary, ...
  'frequency_boundary',est.frequency_boundary,'optimizer_exitflag',est.optimizer_exitflag);
 if ic==1, audit=record; else, audit(ic)=record; end
 if job.preflight
  % Actual production-entry poison test, with fixed observation/calibration data.
  poison=c; poison.Q=999*c.Q; poison.fpeak=999; poison.rho=999*c.rho; poison.vs=999*c.vs;
  pp=uq_public_config(poison,c.nominal_frequency); assert(isequal(public,pp));
  a=p1_compensate(calD,c.offset,public,calD,c.offset,'Offset','CombinedUncertainty');
  b=p1_compensate(calD,c.offset,pp,calD,c.offset,'Offset','CombinedUncertainty');
  assert(isequal(a.direct,b.direct) && isequal(a.split,b.split));
  bad=public; bad.Q=c.Q;
  try
   p1_compensate(calD,c.offset,bad,calD,c.offset,'Offset','CombinedUncertainty');
  catch ex
   rejected_inputs=1; rejected_message=ex.message;
  end
  assert(rejected_inputs==1);
  uq_writejson(fullfile(outdir,'production_entry_poison_test.json'), ...
   struct('actual_entry_called',true,'direct_max_abs_difference',max(abs(a.direct-b.direct),[],'all'), ...
   'split_max_abs_difference',max(abs(a.split-b.split),[],'all'),'extra_truth_field_rejected',true,'message',rejected_message));
 end
 for isu=1:numel(job.suites)
  su=job.suites(isu); est_use=est; cm_use=cm; models_use=invmodels; R_use=R;
  if strcmp(su.kind,'velocity')
   bad=public; bad.vp=public.vp*su.vp_multiplier; bad.h=public.h*su.h_multiplier;
   est_use=p1_calibrate_ensemble(calD,c.offset,bad); cm_use=uq_inverse_config(bad,est_use);
   models_use=prepare_models_v2(c.offset,cm_use); R_use=avo_projector_v2(models_use);
   uq_writejson(fullfile(outdir,sprintf('velocity_%d.json',cseed)), ...
    struct('vp_multiplier',su.vp_multiplier,'h_multiplier',su.h_multiplier, ...
    'original_vp',public.vp,'inverse_vp',bad.vp,'original_h',public.h,'inverse_h',bad.h, ...
    'operator_frobenius_change',norm(models_use{6}.G-invmodels{6}.G,'fro'),'calibration',est_use));
  end
  for isc=1:numel(su.scenarios)
   scname=su.scenarios{isc}; [models,domain,map,scene]=uq_case(scname,c,models_use);
   for seed=reshape(job.test_seeds,1,[])
    [D,trueSigma]=observation(scene,c,trueModels,map,seed,su);
    block=struct([]); nb=0;
    for ib=1:numel(job.methods)
     method=job.methods{ib};
     [out,prior,preset]=p1_solve_branch(D,models,map,cm_use,domain,R_use,method);
     au=uq_audit_outputs(out); flags=flags+au.flags; maxnormal=max(maxnormal,au.max_normal);
     % Numerically identical inverse configuration and data for both routes.
     for ir=1:2
      route='Direct'; S=out.direct; rec=out.direct_record;
      if ir==2, route='Split'; S=out.split; rec=out.split_record; end
      row=uq_metric_row(S,scene,c,c.model_name,su.name,seed,method,route,out,prior,preset,est_use);
      row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma); row.CalibrationSeed=cseed;
      row.Group=job.group; row.InitialFlag=rec.initial_flag; row.Iterations=rec.iterations;
      row.PriorSigmaSum=sum(prior.sigma); row.PriorLambdaSum=sum(prior.lambda);
      row.PriorHFro=norm(prior.H,'fro'); row.CalibrationQBoundary=est_use.q_boundary;
      row.CalibrationFrequencyBoundary=est_use.frequency_boundary;
      row.CalibrationOptimizerFlag=est_use.optimizer_exitflag;
      row.AllSelectionSolverFlags=au.flags; row.AllSelectionMaxNormalResidual=au.max_normal;
      nb=nb+1; if nb==1, block=row; else, block(nb)=row; end
     end
    end
    obs_count=obs_count+1; count=count+nb;
    file=fullfile(outdir,'metrics.csv');
    if obs_count==1, writetable(struct2table(block),file);
    else, writetable(struct2table(block),file,'WriteMode','append','WriteVariableNames',false); end
    fprintf('%s %s %s cal=%d seed=%d observations=%d rows=%d elapsed=%.1fs\n', ...
     c.model_name,su.name,scname,cseed,seed,obs_count,count,toc(start));
    uq_writejson(fullfile(outdir,'progress.json'),struct('observations',obs_count,'rows',count, ...
     'elapsed_seconds',toc(start),'max_normal_residual',maxnormal,'solver_flags',flags));
   end
  end
 end
end
uq_writejson(fullfile(outdir,'calibration_operator_audit.json'),audit);
uq_writejson(fullfile(outdir,'checks.json'),struct('rows',count,'observations',obs_count, ...
 'solver_flags',flags,'max_normal_residual',maxnormal,'elapsed_seconds',toc(start), ...
 'matlab_version',version,'computer',computer,'completed',true));
% Retain all failures; experiment completion and scientific acceptance are separate.
fprintf('COMPLETED rows=%d observations=%d elapsed=%.1fs\n',count,obs_count,toc(start));
end

function [D,sigma]=observation(scene,c,trueModels,map,seed,su)
if any(strcmp(su.kind,{'normal','velocity'}))
 [D,sigma]=uq_observation(scene,c,trueModels,map,seed,su.snr_db); return;
end
[n,nx]=size(scene.truth); clean=zeros(n,nx); N=c.nfft;
for j=1:nx
 if strcmp(su.kind,'phase')
  rad=su.phase_deg*pi/180;
  phase=[1;exp(-1i*rad)*ones(N/2-1,1);1;exp(1i*rad)*ones(N/2-1,1)];
  wr=real(ifft(fft(scene.truth(:,j),N).*phase)); clean(:,j)=trueModels{map(j)}.G*wr(1:n);
 else
  raw=trueModels{map(j)}.G*scene.truth(:,j);
  if strcmp(su.kind,'nmo')
   f=(0:N/2)'/(N*c.dt); wf=fft(raw,N); shift=su.max_nmo_shift_s*((j-1)/(nx-1))^2;
   half=wf(1:N/2+1).*exp(-2i*pi*f*shift);
   padded=real(ifft([half;conj(half(end-1:-1:2))])); raw=padded(1:n);
  end
  clean(:,j)=raw;
 end
end
coherent=zeros(n,nx);
if strcmp(scene.name,'CoherentNoiseAVO')
 for j=1:nx, z=pi*35*(c.t-(mean(c.event_t(1:2))+.004*(j-1))); coherent(:,j)=(1-2*z.^2).*exp(-z.^2); end
 coherent=.25*norm(clean,'fro')*coherent/norm(coherent,'fro');
end
rng(seed,'twister'); noise=randn(n,nx);
if strcmp(su.kind,'colored')
 rho=su.ar1_rho;
 for t=2:n, noise(t,:)=rho*noise(t-1,:)+sqrt(1-rho^2)*noise(t,:); end
end
noise=noise./sqrt(mean(noise.^2,1)); sigma=sqrt(mean(clean.^2,1))*10^(-su.snr_db/20);
D=clean+noise.*sigma+coherent;
end
function restore(p,r), path(p); rng(r); end
