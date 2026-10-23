function run_universal_screen(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));start=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','universal_plan.json')));
uq_writejson(fullfile(outdir,'executed_plan.json'),plan);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
cache=plan.cache;files=dir(fullfile(cache,'cases','*.mat'));
cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
rows=struct([]);count=0;maxnormal=0;replay=0;previous='';
for k=1:numel(files)
 z=load(fullfile(files(k).folder,files(k).name));meta=z.meta;
 if ~strcmp(previous,meta.model)
  a=jsondecode(fileread(fullfile(cache,['calibration_' meta.model '.json'])));
  c=config_model();c.model_name=meta.model;c.nominal_frequency=60;
  if ~strcmp(meta.model,'OriginalReference')
   list=cached.additional_models;ix=find(strcmp({list.name},meta.model),1);m=list(ix);
   for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end
   c.h=diff([0,c.event_t]).*c.vp(1:3)/2;c.nominal_frequency=m.nominal_frequency;
  end
  previous=meta.model;previousDomain='';
 end
 if ~strcmp(previousDomain,meta.domain)
  est=p1_candidate_domain_estimate(a.diagnostics.fixed_estimate,a.estimate,a.corrected,policy,meta.domain);
  cm=uq_inverse_config(a.corrected,est);cm.candidate_profile=policy.selected_profile;cm.joint_prepare_only=true;
  if strcmp(meta.domain,'Offset'),map=1:11;models=prepare_models_v2(z.offsets,cm);R=avo_projector_v2(models);
  else,map=6*ones(1,11);one=prepare_models_v2(z.offsets(1),cm);models=repmat(one,1,11);R=sparse(c.nt*11,c.nt*11);end
  previousDomain=meta.domain;
 end
 D=z.D;[~,~,~,ctx]=p1_solve_branch(D,models,map,cm,meta.domain,R,'CombinedUncertainty');
 scene=make_stability_scene(meta.scenario,c);assert(norm(scene.truth-z.truth,'fro')<1e-12);
 for im=1:3
  label='Baseline';S=z.Snew;normal=meta.normal_new;
  if im>1
   label=plan.candidates{im-1};out=universal_candidate(ctx,cm,label);S=out.split;
   audit=uq_audit_outputs(out);assert(audit.flags==0&&audit.max_normal<2e-7);normal=audit.max_normal;
   maxnormal=max(maxnormal,normal);
  end
  score=score_stability(S,scene,c);
  if im==1,replay=max(replay,abs(score.RelativeError-meta.new_score.RelativeError));end
  row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method',label);
  f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end
  row.NormalResidual=normal;count=count+1;if count==1,rows=row;else,rows(count)=row;end
  if im==2,mixture=S;elseif im==3,sourceonly=S;end
 end
 baseline=z.Snew;truth=z.truth;t=z.t;observed=D;
 save(fullfile(outdir,'cases',files(k).name),'mixture','sourceonly','baseline','truth','t','observed','meta','-v7');
 writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
 uq_writejson(fullfile(outdir,'progress.json'),struct('completed_cases',k,'total',numel(files),'elapsed_seconds',toc(start)));
 fprintf('%d/%d %s elapsed=%.1fs\n',k,numel(files),meta.id,toc(start));
end
assert(numel(files)==60&&count==180&&replay<1e-12);
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',60,'rows',180,'max_normal',maxnormal,'baseline_score_replay',replay,'elapsed_seconds',toc(start)));
end
