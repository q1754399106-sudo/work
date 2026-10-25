function run_calibration_shrinkage_diagnostic(outdir)
% Development-only shrinkage of uncertain calibration parameters.
% Every candidate is still a complete inverse solve with its own physical G.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);start=tic;
pointer=jsondecode(fileread(fullfile(root,'inputs','section_diagnostic_source.json')));source=pointer.run;
plan=jsondecode(fileread(fullfile(source,'matlab','executed_plan.json')));m=plan.additional_models(1);c=config_model();
for key={'vp','vs','rho','Q','event_t','fpeak'},k=key{1};c.(k)=reshape(m.(k),1,[]);end
c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.model_name=m.name;c.nominal_frequency=m.nominal_frequency;
public=uq_public_config(c,c.nominal_frequency);cseed=plan.calibration_seeds(1);seed=plan.test_seeds(1);
file=fullfile(source,'matlab',sprintf('calibration_%s_Main_%d.json',c.model_name,cseed));d=jsondecode(fileread(file));
fixed=d.fixed_estimate;layered=d.layered_estimate;corrected=public;corrected.vp=public.vp*d.actual_velocity_scale;
scene=make_stability_scene('WeakSection',c);map=6*ones(1,11);R=sparse(c.nt*11,c.nt*11);
[aa,bb]=ray_absorption(c.t,c.offset(6),c);trueModels=repmat({struct('G',attenuation_matrix(aa+bb,c))},1,11);
weights=[0,.25,.5,1];profiles=struct('name','OriginalPrior','confidence_floor',1,'section_modal_blend',1,'section_target',-1,'section_gamma',10);
profiles(2)=struct('name','LightPrior','confidence_floor',.2,'section_modal_blend',.25,'section_target',1,'section_gamma',3);
uq_writejson(fullfile(outdir,'executed_plan.json'),struct('development_only',true,'model',m,'source_calibration',file,'weights',weights,'profiles',profiles,'calibration_seed',cseed,'seed',seed));
uq_writejson(fullfile(outdir,'reused_calibration.json'),d);rows=struct([]);count=0;flags=0;maxnormal=0;
for weight=weights
 est=fixed;
 if weight==1,est=layered;
 elseif weight>0
  est.layer_Q=exp((1-weight)*log(fixed.effective_Q)+weight*log(reshape(layered.layer_Q,1,[])));
  est.estimated_frequency=exp((1-weight)*log(fixed.estimated_frequency)+weight*log(layered.estimated_frequency));
  travel=2*corrected.h./corrected.vp(1:3);est.effective_Q=sum(travel)/sum(travel./est.layer_Q(1:3));
 end
 cm=uq_inverse_config(corrected,est);one=prepare_models_v2(c.offset(6),cm);mm=repmat(one,1,11);
 for isnr=1:2
  snr=18;suite='Main';if isnr==2,snr=12;suite='LowerSNR';end
  [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,snr);
  for ip=1:numel(profiles)
   local=cm;local.candidate_profile=profiles(ip);method=sprintf('Shrink%.2f_%s',weight,profiles(ip).name);
   [out,prior,preset]=p1_solve_branch(D,mm,map,local,'Section',R,'CombinedUncertainty');quality=uq_audit_outputs(out);
   flags=flags+quality.flags;maxnormal=max(maxnormal,quality.max_normal);
   for ir=1:2
    route='Direct';S=out.direct;if ir==2,route='Split';S=out.split;end
    row=uq_metric_row(S,scene,c,c.model_name,suite,seed,method,route,out,prior,preset,est);
    row.CalibrationSeed=cseed;row.CalibrationBlend=weight;row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma);row.AllSelectionSolverFlags=quality.flags;row.AllSelectionMaxNormalResidual=quality.max_normal;
    count=count+1;if count==1,rows=row;else,rows(count)=row;end
   end
   uq_writejson(fullfile(outdir,[suite '_' method '_waveforms.json']),struct('t',c.t,'observed',D,'truth',scene.truth,'direct',out.direct,'split',out.split));
   writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
   fprintf('%s %s rows=%d elapsed=%.1fs normal=%.3e flags=%d\n',suite,method,count,toc(start),maxnormal,flags);
  end
 end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'rows',count,'solver_flags',flags,'max_normal_residual',maxnormal,'development_only',true));
end
function restore(p,r),path(p);rng(r);end
