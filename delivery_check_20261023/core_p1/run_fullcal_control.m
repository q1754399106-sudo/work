function run_fullcal_control(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));cache=plan.cache;
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
assert(~isfolder(outdir));mkdir(outdir);start=tic;files=dir(fullfile(cache,'cases','*.mat'));rows=struct([]);previous='';maxnormal=0;
for k=1:numel(files)
 z=load(fullfile(files(k).folder,files(k).name));meta=z.meta;
 if ~strcmp(previous,meta.model)
  a=jsondecode(fileread(fullfile(cache,['calibration_' meta.model '.json'])));for key={'vp','h','offset','lambda_grid'},f=key{1};a.corrected.(f)=reshape(a.corrected.(f),1,[]);end
  a.corrected.t=a.corrected.t(:);cm=uq_inverse_config(a.corrected,a.estimate);cm.candidate_profile=policy.selected_profile;cm.joint_prepare_only=true;
  one=prepare_models_v2(a.corrected.offset(6),cm);models=repmat(one,1,11);R=sparse(cm.nt*11,cm.nt*11);map=6*ones(1,11);previous=meta.model;
  c=config_model();c.model_name=meta.model;if ~strcmp(meta.model,'OriginalReference'),list=cached.additional_models;ix=find(strcmp({list.name},meta.model),1);m=list(ix);for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end;c.h=diff([0,c.event_t]).*c.vp(1:3)/2;end
 end
 S=z.Snew;normal=meta.normal_new;
 if strcmp(meta.domain,'Section')
  [~,~,~,ctx]=p1_solve_branch(z.D,models,map,cm,'Section',R,'CombinedUncertainty');out=joint_select(ctx,'DecoupledWeak',.3,{});S=out.split;au=uq_audit_outputs(out);normal=au.max_normal;assert(au.flags==0);
 end
 scene=make_stability_scene(meta.scenario,c);score=score_stability(S,scene,c);row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method','FullCalLinear');
 f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end;row.NormalResidual=normal;
 if k==1,rows=row;else,rows(k)=row;end;maxnormal=max(maxnormal,normal);writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
 fprintf('%d/60 %s error=%.5f elapsed=%.1fs\n',k,meta.id,score.RelativeError,toc(start));
end
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',60,'max_normal',maxnormal,'role','Post-selection attribution control, not used to select the frozen candidate','elapsed_seconds',toc(start)));
end
