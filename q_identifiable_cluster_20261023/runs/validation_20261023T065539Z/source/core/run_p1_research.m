function run_p1_research(action,outdir)
% Shared-session P1 research entry. Never touches base workspace, path, or RNG permanently.
root=fileparts(mfilename('fullpath')); pkg=fileparts(root);
assert(~isfolder(outdir),'Output must be a new directory');
oldpath=path; oldrng=rng; cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));
switch action
 case {'validate','smoke'}, validate_p1_unknown_q(outdir);
 case {'regression','replay'}, run_research(action,outdir);
 case {'p1_experiment','full'}, run_p1_grid(outdir,pkg);
 otherwise, error('Unknown action: %s',action);
end
end
function restore(p,r), path(p); rng(r); end
function run_p1_grid(outdir,pkg)
mkdir(outdir); start=tic;
spec=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'inputs','p1_experiment_spec.json')));
splits={'historical','development','heldout_new_geology'};
rows_hist=struct([]); n_hist=0;
rows_all=struct([]); n_all=0;
rows_stress=struct([]); n_stress=0;
flags=0; maxnormal=0;
for is=1:numel(splits)
 sname=splits{is}; [configs,~]=p1_models(sname);
 for im=1:numel(configs)
  c=configs{im}; public=uq_public_config(c,c.nominal_frequency);
  trueModels=prepare_models_v2(c.offset,c);
  [~,~,~,calScene]=uq_case('SmoothAVO',c,trueModels);
  for ic=1:numel(spec.calibration_seeds)
   cseed=spec.calibration_seeds(ic);
   [calD,~]=uq_observation(calScene,c,trueModels,1:11,cseed,18);
   est=p1_calibrate_ensemble(calD,c.offset,public);
   cm=uq_inverse_config(public,est);
   invmodels=prepare_models_v2(c.offset,cm);
   R_avo=avo_projector_v2(invmodels);
   if strcmp(sname,'historical') && cseed==20261030
    for isu=1:2
     suite='Main'; snr=18; scenarios=spec.scenarios;
     if isu==2, suite='LowerSNR'; snr=12; scenarios=spec.stress_scenarios; end
     for isc=1:numel(scenarios)
      scname=char(scenarios{isc});
      [models,domain,map,scene]=uq_case(scname,c,invmodels);
      for seed=reshape(spec.test_seeds,1,[])
       [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,snr);
       for ib=1:numel(spec.methods)
        bname=spec.methods(ib).name;
        [out,prior,preset]=p1_solve_branch(D,models,map,cm,domain,R_avo,bname);
        au=uq_audit_outputs(out); flags=flags+au.flags; maxnormal=max(maxnormal,au.max_normal);
        for ir=1:2
         route='Direct'; S=out.direct; if ir==2, route='Split'; S=out.split; end
         row=uq_metric_row(S,scene,c,c.model_name,suite,seed,bname,route,out,prior,preset,est);
         row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma); row.CalibrationSeed=cseed;
         n_hist=n_hist+1; if n_hist==1, rows_hist=row; else, rows_hist(n_hist)=row; end
         n_all=n_all+1; if n_all==1, rows_all=row; else, rows_all(n_all)=row; end
        end
       end
      end
     end
    end
    % Stress suites on historical models (SmoothAVO & CoherentNoiseAVO)
    st_names={'StressNonRickerPhase','StressVelocityError','StressResidualNMO','StressColoredNoise'};
    st_scs={'SmoothAVO','CoherentNoiseAVO'};
    for ist=1:numel(st_names)
     stname=st_names{ist};
     cm_st=cm; invmodels_st=invmodels; R_avo_st=R_avo; est_st=est;
     if strcmp(stname,'StressVelocityError')
      public_st=public; public_st.vp=1.04*public.vp; public_st.h=1.035*public.h;
      est_st=p1_calibrate_ensemble(calD,c.offset,public_st);
      cm_st=uq_inverse_config(public_st,est_st);
      invmodels_st=prepare_models_v2(c.offset,cm_st);
      R_avo_st=avo_projector_v2(invmodels_st);
     end
     for isc=1:numel(st_scs)
      scname=st_scs{isc};
      [models_st,domain_st,map_st,scene_st]=uq_case(scname,c,invmodels_st);
      [D_st,trueSigma_st]=synthesize_stress_obs(scene_st,c,trueModels,map_st,20261102,18,stname);
      for ib=1:numel(spec.methods)
       bname=spec.methods(ib).name;
       [out,prior,preset]=p1_solve_branch(D_st,models_st,map_st,cm_st,domain_st,R_avo_st,bname);
       au=uq_audit_outputs(out); flags=flags+au.flags; maxnormal=max(maxnormal,au.max_normal);
       for ir=1:2
        route='Direct'; S=out.direct; if ir==2, route='Split'; S=out.split; end
        row=uq_metric_row(S,scene_st,c,c.model_name,stname,20261102,bname,route,out,prior,preset,est_st);
        row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma_st); row.CalibrationSeed=cseed;
        n_stress=n_stress+1; if n_stress==1, rows_stress=row; else, rows_stress(n_stress)=row; end
       end
      end
     end
    end
   else
    if strcmp(sname,'historical')
     suite_mc='MultiCalMain'; seed_mc=20261101+mod(cseed,5);
    elseif strcmp(sname,'development')
     suite_mc='DevMain'; seed_mc=20261101+floor((cseed-20261030)/5);
    else
     suite_mc='NewGeologyMain'; seed_mc=20261101+floor((cseed-20261030)/5);
    end
    for isc=1:numel(spec.scenarios)
     scname=char(spec.scenarios{isc});
     [models,domain,map,scene]=uq_case(scname,c,invmodels);
     [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed_mc,18);
     for ib=1:numel(spec.methods)
      bname=spec.methods(ib).name;
      [out,prior,preset]=p1_solve_branch(D,models,map,cm,domain,R_avo,bname);
      au=uq_audit_outputs(out); flags=flags+au.flags; maxnormal=max(maxnormal,au.max_normal);
      for ir=1:2
       route='Direct'; S=out.direct; if ir==2, route='Split'; S=out.split; end
       row=uq_metric_row(S,scene,c,c.model_name,suite_mc,seed_mc,bname,route,out,prior,preset,est);
       row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma); row.CalibrationSeed=cseed;
       n_all=n_all+1; if n_all==1, rows_all=row; else, rows_all(n_all)=row; end
      end
     end
    end
   end
  end
 end
end
writetable(struct2table(rows_hist),fullfile(outdir,'metrics.csv'));
writetable(struct2table(rows_hist),fullfile(outdir,'metrics_p1_historical_paired.csv'));
writetable(struct2table(rows_all),fullfile(outdir,'metrics_p1_all_models_multical.csv'));
writetable(struct2table(rows_stress),fullfile(outdir,'metrics_p1_stress_suites.csv'));
uq_writejson(fullfile(outdir,'checks.json'),struct('solver_flags',flags,'max_normal_residual',maxnormal,'rows_historical',n_hist,'rows_all',n_all,'rows_stress',n_stress,'elapsed_seconds',toc(start)));
end
function [D,trueSigma]=synthesize_stress_obs(scene,c,trueModels,map,seed,snr,stname)
[n,nx]=size(scene.truth); clean=zeros(n,nx);
if strcmp(stname,'StressNonRickerPhase')
 N=c.nfft; rad=30*pi/180;
 phase=[1; exp(-1i*rad)*ones(N/2-1,1); 1; exp(1i*rad)*ones(N/2-1,1)];
 for j=1:nx
  wf=fft(scene.truth(:,j),N); w_rot=real(ifft(wf.*phase,N));
  clean(:,j)=trueModels{map(j)}.G*w_rot(1:n);
 end
elseif strcmp(stname,'StressResidualNMO')
 N=c.nfft; f=(0:N/2)'/(N*c.dt);
 for j=1:nx
  raw_j=trueModels{map(j)}.G*scene.truth(:,j);
  dt_shift=0.004*((j-1)/10)^2;
  wf=fft(raw_j,N); half=wf(1:N/2+1).*exp(-2i*pi*f*dt_shift);
  full_sp=[half; conj(half(end-1:-1:2))];
  clean(:,j)=real(ifft(full_sp,N)); clean(:,j)=clean(1:n,j);
 end
else
 for j=1:nx, clean(:,j)=trueModels{map(j)}.G*scene.truth(:,j); end
end
coherent=zeros(n,nx);
if strcmp(scene.name,'CoherentNoiseAVO')
 for j=1:nx
  z=pi*35*(c.t-(mean(c.event_t(1:2))+0.004*(j-1)));
  coherent(:,j)=(1-2*z.^2).*exp(-z.^2);
 end
 coherent=0.25*norm(clean,'fro')*coherent/norm(coherent,'fro');
end
rng(seed,'twister'); raw=randn(n,nx);
if strcmp(stname,'StressColoredNoise')
 rho=0.55; col=zeros(n,nx); col(1,:)=raw(1,:); sf=sqrt(1-rho^2);
 for t=2:n, col(t,:)=rho*col(t-1,:)+sf*raw(t,:); end
 raw=col;
end
noise=raw./sqrt(mean(raw.^2,1));
trueSigma=sqrt(mean(clean.^2,1))*10^(-snr/20);
D=clean+noise.*trueSigma+coherent;
end
