function replay_confirmation_case(runroot,outdir)
% Independent production-API replay of a predetermined prospective case.
% No changes to the frozen numerical method or acceptance are made here.
oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
core=fullfile(runroot,'source','core_p1');addpath(core);addpath(fullfile(core,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);plan=jsondecode(fileread(fullfile(runroot,'FROZEN_EXECUTION_PLAN.json')));
index=find(~[plan.jobs.reuse_original],1,'first');job=plan.jobs(index);m=job.model;
c=config_model();for key={'vp','vs','rho','Q','event_t','fpeak'},k=key{1};c.(k)=reshape(m.(k),1,[]);end
c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;
cseed=job.calibration_seeds(1);seed=job.test_seeds(1);
saved=load(fullfile(runroot,m.name,'matlab',sprintf('calibration_input_%d.mat',cseed)),'calD','public');
public=saved.public;calD=saved.calD;scene=make_stability_scene('WeakSection',c);map=6*ones(1,11);
[a,b]=ray_absorption(c.t,c.offset(6),c);trueModels=repmat({struct('G',attenuation_matrix(a+b,c))},1,11);
isu=find(strcmp({job.suites.name},'Main'),1);su=job.suites(isu);
[D,trueSigma]=p1_validation_observation(scene,c,trueModels,map,seed,su);
out=combined_uncertainty_v2(D,scene.offsets,public,calD,c.offset,'Section');
quality=uq_audit_outputs(out);est=out.estimated_Q;prior=out.prior;preset=out.preset;preset.name='CombinedUncertaintyV2';
corrected=public;corrected.vp=public.vp*out.estimated_velocity_scale;cm=uq_inverse_config(corrected,est);
rows=struct([]);
for ir=1:2
 route='Direct';S=out.direct;rec=out.direct_record;if ir==2,route='Split';S=out.split;rec=out.split_record;end
 row=uq_metric_row(S,scene,c,c.model_name,'Main',seed,'CombinedUncertaintyV2',route,out,prior,preset,est);
 row.CalibrationSeed=cseed;row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma);row.Group=job.group;
 row.InitialFlag=rec.initial_flag;row.Iterations=rec.iterations;row.PriorSigmaSum=sum(prior.sigma);row.PriorLambdaSum=sum(prior.lambda);row.PriorHFro=norm(prior.H,'fro');
 row.CalibrationQBoundary=est.q_boundary;row.CalibrationFrequencyBoundary=est.frequency_boundary;row.CalibrationOptimizerFlag=est.optimizer_exitflag;
 row.AllSelectionSolverFlags=quality.flags;row.AllSelectionMaxNormalResidual=quality.max_normal;
 row.VelocityCorrection=mean(cm.vp./public.vp);row.LayerQ1=cm.Q(1);row.LayerQ2=cm.Q(2);row.LayerQ3=cm.Q(3);row.LayerQ4=cm.Q(4);
 row.ComputationAliasOf='';row.CompensationApplied=1;row.DataEnergy=sum(D(:).^2);row.RowProvenance='new_MATLAB_execution';
 if ir==1,rows=row;else,rows(ir)=row;end
end
writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
example=jsondecode(fileread(fullfile(runroot,m.name,'matlab','example_WeakSection.json')));
assert(isequal(out.split(:,6),example.CombinedUncertaintyV2(:)),'Full precision serialized waveform differs from original');
assert(quality.flags==0&&quality.max_normal<=2e-7);
uq_writejson(fullfile(outdir,'checks.json'),struct('passed',true,'model',m.name,'calibration_seed',cseed,'observation_seed',seed, ...
 'scenario','WeakSection','suite','Main','routes',{{'Direct','Split'}},'production_api_called',true,'fixed_trace6_waveform_exact',true,'selection_rule','First prospective model, first calibration and observation seed, Main WeakSection predetermined before inspecting confirmation outcomes'));
fprintf('PASS production API replay, both routes and exact trace waveform\n');
end
function restore(p,r),path(p);rng(r);end
