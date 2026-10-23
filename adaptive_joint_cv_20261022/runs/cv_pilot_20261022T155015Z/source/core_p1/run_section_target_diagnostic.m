function run_section_target_diagnostic(outdir,mode)
% Development diagnosis: separate discrepancy strength from structure damping.
% Reuse a completed calibration, never calibrate against reconstruction truth.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
if nargin<2,mode='target';end
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
pointer=jsondecode(fileread(fullfile(root,'inputs','section_diagnostic_source.json')));
prior=pointer.run;plan=jsondecode(fileread(fullfile(prior,'matlab','executed_plan.json')));
m=plan.additional_models(1);c=config_model();
for key={'vp','vs','rho','Q','event_t','fpeak'},k=key{1};c.(k)=reshape(m.(k),1,[]);end
c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;
public=uq_public_config(c,c.nominal_frequency);cseed=plan.calibration_seeds(1);seed=plan.test_seeds(1);
file=fullfile(prior,'matlab',sprintf('calibration_%s_Main_%d.json',c.model_name,cseed));d=jsondecode(fileread(file));
est=d.layered_estimate;corrected=public;corrected.vp=public.vp*d.actual_velocity_scale;
cm=uq_inverse_config(corrected,est);models=prepare_models_v2(c.offset,cm);R=avo_projector_v2(models);
trueModels=cell(1,11);for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);trueModels{j}=struct('G',attenuation_matrix(a+b,c));end
profiles=struct('name','ReplayConfidenceIdentityPenalty','confidence_floor',.2,'section_modal_blend',0,'section_target',1,'section_gamma',3);
profiles(2)=struct('name','Target085Gamma3','confidence_floor',.2,'section_modal_blend',0,'section_target',.85,'section_gamma',3);
profiles(3)=struct('name','Target075Gamma3','confidence_floor',.2,'section_modal_blend',0,'section_target',.75,'section_gamma',3);
profiles(4)=struct('name','Target085Gamma0','confidence_floor',1,'section_modal_blend',0,'section_target',.85,'section_gamma',0);
profiles(5)=struct('name','Target100Gamma0','confidence_floor',1,'section_modal_blend',0,'section_target',1,'section_gamma',0);
if strcmp(mode,'sourceband')
 profiles=repmat(struct('name','','confidence_floor',1,'section_modal_blend',0,'section_target',1,'section_gamma',0,'sourceband_blend',1,'sourceband_floor',.05),1,4);
 profiles(1).name='SourceBand005Gamma0';profiles(2).name='SourceBand020Gamma0';profiles(2).sourceband_floor=.2;
 profiles(3).name='SourceBand005Gamma3';profiles(3).section_gamma=3;profiles(3).confidence_floor=.2;
 profiles(4).name='SourceBand020Gamma3';profiles(4).sourceband_floor=.2;profiles(4).section_gamma=3;profiles(4).confidence_floor=.2;
end
uq_writejson(fullfile(outdir,'executed_plan.json'),struct('development_only',true,'source_calibration',file,'profiles',profiles,'model',m,'calibration_seed',cseed,'seed',seed,'scenario','WeakSection'));
rows=struct([]);count=0;flags=0;maxnormal=0;
[mm,domain,map,scene]=uq_case('WeakSection',c,models);
for isnr=1:2
 snr=18;suite='Main';if isnr==2,snr=12;suite='LowerSNR';end
 [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,snr);
 for ip=1:numel(profiles)
  local=cm;local.candidate_profile=profiles(ip);
  [out,priorout,preset]=p1_solve_branch(D,mm,map,local,domain,R,'CombinedUncertainty');
  quality=uq_audit_outputs(out);flags=flags+quality.flags;maxnormal=max(maxnormal,quality.max_normal);
  for ir=1:2
   route='Direct';S=out.direct;if ir==2,route='Split';S=out.split;end
   row=uq_metric_row(S,scene,c,c.model_name,suite,seed,profiles(ip).name,route,out,priorout,preset,est);
   row.CalibrationSeed=cseed;row.EstimatedSigmaRatio=mean(priorout.sigma./trueSigma);row.AllSelectionSolverFlags=quality.flags;row.AllSelectionMaxNormalResidual=quality.max_normal;
   count=count+1;if count==1,rows=row;else,rows(count)=row;end
  end
  waveform=struct('t',c.t,'observed',D,'truth',scene.truth,'direct',out.direct,'split',out.split,'sigma',priorout.sigma,'scale',out.scale,'data_residual',out.direct_record.data_residual,'target',out.target);
  uq_writejson(fullfile(outdir,[suite '_' profiles(ip).name '_waveforms.json']),waveform);
  writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
  fprintf('%s %s rows=%d elapsed=%.1fs maxnormal=%.3e flags=%d\n',suite,profiles(ip).name,count,toc(start),maxnormal,flags);
 end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'rows',count,'solver_flags',flags,'max_normal_residual',maxnormal,'development_only',true));
end
function restore(p,r),path(p);rng(r);end
