function run_universal_consensus(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
assert(~isfolder(outdir));mkdir(outdir);mkdir(fullfile(outdir,'cases'));
plan=jsondecode(fileread(fullfile(root,'inputs','consensus_run.json')));
previous='';rows=struct([]);files=dir(fullfile(plan.source,'cases','*.mat'));
cached=jsondecode(fileread(fullfile(root,'inputs','cached_shrinkage','trial_plan.json')));
for k=1:numel(files)
 z=load(fullfile(files(k).folder,files(k).name));meta=z.meta;
 if ~strcmp(previous,meta.model)
  c=config_model();c.model_name=meta.model;
  if ~strcmp(meta.model,'OriginalReference')
   list=cached.additional_models;ix=find(strcmp({list.name},meta.model),1);m=list(ix);
   for key={'vp','vs','rho','Q','event_t','fpeak'},f=key{1};c.(f)=reshape(m.(f),1,[]);end
   c.h=diff([0,c.event_t]).*c.vp(1:3)/2;
  end
  previous=meta.model;
 end
 scene=make_stability_scene(meta.scenario,c);assert(norm(scene.truth-z.truth,'fro')<1e-12);
 S=.5*z.mixture+.5*z.sourceonly;score=score_stability(S,scene,c);
 row=struct('CaseID',meta.id,'Model',meta.model,'Domain',meta.domain,'Scenario',meta.scenario,'SNR_dB',meta.snr_db,'Seed',meta.seed,'Method','Consensus');
 f=fieldnames(score);for i=1:numel(f),row.(f{i})=score.(f{i});end
 row.NormalResidual=NaN; % No single normal equation for an estimator ensemble.
 if k==1,rows=row;else,rows(k)=row;end
 t=z.t;truth=z.truth;baseline=z.baseline;observed=z.observed;
 save(fullfile(outdir,'cases',files(k).name),'S','t','truth','baseline','observed','meta','-v7');
end
assert(numel(files)==60);writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'cases',60,'weights',[.5,.5],'source',plan.source));
end
