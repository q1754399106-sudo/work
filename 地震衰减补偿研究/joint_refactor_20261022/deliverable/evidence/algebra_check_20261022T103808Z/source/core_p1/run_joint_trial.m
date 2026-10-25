function run_joint_trial(outdir)
% Paired development experiment; never selects methods using synthesis truth.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','joint_plan.json')));uq_writejson(fullfile(outdir,'executed_plan.json'),plan);
[configs,~]=p1_models('development');cached=fullfile(root,'inputs','cached_shrinkage');
extra=jsondecode(fileread(fullfile(cached,'trial_plan.json')));
for im=1:numel(extra.additional_models)
 m=extra.additional_models(im);c=config_model();for f={'vp','vs','rho','Q','event_t','fpeak'},k=f{1};c.(k)=reshape(m.(k),1,[]);end
 c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;configs{end+1}=c;
end
rows=struct([]);count=0;maxnormal=0;flags=0;identity=0;replay=NaN;shared_difference=0;
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
for im=1:numel(configs)
 c=configs{im};public=uq_public_config(c,c.nominal_frequency);
 cal=jsondecode(fileread(fullfile(cached,['calibration_' c.model_name '.json'])));
 corrected=public;corrected.vp=public.vp*cal.actual_velocity_scale;
 trueModels=cell(1,11);for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(a+b,c));end
 for ic=1:numel(plan.calibrations)
  calibration=plan.calibrations{ic};
  for id=1:2
   domain='Offset';if id==2,domain='Section';end
   if strcmp(calibration,'Old'),est=cal.fixed_estimate;provided=public;
   else,est=p1_candidate_domain_estimate(cal.fixed_estimate,cal.layered_estimate,corrected,policy,domain);provided=corrected;end
   cm=uq_inverse_config(provided,est);allmodels=prepare_models_v2(c.offset,cm);Ravo=avo_projector_v2(allmodels);
   offsets=c.offset;if id==2,offsets=repmat(c.offset(6),1,11);end
   ensemble=joint_physical_ensemble(cm,offsets);
   for j=1:11,identity=max(identity,ensemble{j}.identity_error);end
   % Same profile for old/new calibration: calibration is the only crossed factor.
   cm.candidate_profile=policy.selected_profile;
   for isc=1:numel(plan.scenarios)
    name=plan.scenarios{isc};if endsWith(name,'Section')~=(id==2),continue;end
    [models,~,map,scene]=uq_case(name,c,allmodels);
    for seed=reshape(plan.test_seeds,1,[])
     [D,~]=uq_observation(scene,c,trueModels,map,seed,plan.snr_db);
     cm.joint_prepare_only=true;[~,prior,preset,ctx]=p1_solve_branch(D,models,map,cm,domain,Ravo,'CombinedUncertainty');
     for ik=1:numel(plan.methods)
      method=plan.methods{ik};factor=plan.beta_factors(ik);
      if strcmp(method,'Legacy')
       cm.joint_prepare_only=false;[out,~,~]=p1_solve_branch(D,models,map,cm,domain,Ravo,'CombinedUncertainty');
       if isnan(replay)
        test=joint_select(ctx,'CoupledReplay',1,ensemble);
        replay=max(norm(test.direct-out.direct,'fro')/norm(out.direct,'fro'),norm(test.split-out.split,'fro')/norm(out.split,'fro'));
        assert(replay<1e-4,'Refactored solver replay mismatch');
        uq_writejson(fullfile(outdir,'replay_check.json'),struct('relative_difference',replay));
       end
      else,out=joint_select(ctx,method,factor,ensemble);end
      a=uq_audit_outputs(out);flags=flags+a.flags;maxnormal=max(maxnormal,a.max_normal);
      for ir=1:2
       route='Direct';S=out.direct;rec=out.direct_record;if ir==2,route='Split';S=out.split;rec=out.split_record;end
       row=uq_metric_row(S,scene,c,c.model_name,'Main',seed,method,route,out,prior,preset,est);
       row.Calibration=calibration;row.CalibrationSeed=plan.calibration_seed;row.BetaFactor=factor;
       row.AllSelectionFlags=a.flags;row.AllSelectionMaxNormal=a.max_normal;
       row.HRms=norm(ctx.H,'fro')/sqrt(max(1,sum(full(sum(ctx.H.^2,2))>0)));
       row.ModalPenalty=NaN;row.StructurePenalty=NaN;row.QPenalty=NaN;
       if isfield(rec,'modal_penalty'),row.ModalPenalty=rec.modal_penalty;row.StructurePenalty=rec.structure_penalty;row.QPenalty=rec.q_penalty;end
       row.SharedRegularizerDifference=NaN;if any(strcmp(method,{'SharedIsotropic','Physical'})),row.SharedRegularizerDifference=out.regularizer_difference;assert(out.regularizer_difference<1e-14);end
       count=count+1;if count==1,rows=row;else,rows(count)=row;end
      end
     end
     writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
     fprintf('%s %s %s seed=%d rows=%d elapsed=%.1fs normal=%.2e\n',c.model_name,calibration,name,seed,count,toc(start),maxnormal);
    end
   end
  end
 end
end
assert(flags==0);assert(maxnormal<2*structure_options_tolerance());
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'rows',count,'solver_flags',flags, ...
 'max_normal_residual',maxnormal,'expectation_identity_error',identity,'legacy_replay_difference',replay, ...
 'elapsed_seconds',toc(start),'development_only',true,'matlab_version',version));
end
function tol=structure_options_tolerance(),o=structure_options();tol=o.pcg_tolerance;end
function restore(p,r),path(p);rng(r);end
