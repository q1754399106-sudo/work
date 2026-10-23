function run_unknown_q_test(outdir)
assert(~isfolder(outdir)); mkdir(outdir); start=tic; root=fileparts(fileparts(mfilename('fullpath')));
[configs,~]=stability_models(); o=uq_spec(); selected=jsondecode(fileread(fullfile(root,'results','matlab_search2','selected_presets.json')));
rows=struct([]); evidence=struct([]); count=0; inputs=0; flags=0; initialflags=0; fallbacks=0; maxres=0; maxfit=0; endpoints=0;
uq_writejson(fullfile(outdir,'selected_before_test.json'),selected); uq_writejson(fullfile(outdir,'spec_used.json'),o);
for im=2:numel(configs)
 c=configs{im}; p=uq_public_config(c,o.nominal_frequency_hz(im-1)); est=jsondecode(fileread(fullfile(root,'results','matlab_calibrate',['estimate_' c.model_name '.json'])));
 cm=uq_inverse_config(p,est); allmodels=prepare_models_v2(c.offset,cm); trueModels=prepare_models_v2(c.offset,c);
 ix=find(strcmp({selected.Model},c.model_name)); presets=[o.presets(1);o.presets(2);selected(ix).Preset]; methods={'Baseline','Frozen','Tuned'};
 for isu=1:2
  suite='Main'; snr=18; scenarios=o.test_scenarios;
  if isu==2, suite='LowerSNR'; snr=12; scenarios=o.stress_scenarios; end
  for ic=1:numel(scenarios)
   name=char(scenarios{ic}); [models,domain,map,scene]=uq_case(name,c,allmodels);
   for seed=reshape(o.test_seeds,1,[])
    [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,snr); Z=zeros(c.nt,11,2,3);
    for ip=1:3
     preset=presets(ip);
     same=find(arrayfun(@(pp)isequal(pp,preset),presets(1:ip-1)),1);
     if ~isempty(same), out=outs{same}; prior=priors{same};
     else, prior=uq_build_prior(D,models,cm,domain,preset); out=uq_select(D,models,prior,preset); end
     outs{ip}=out; priors{ip}=prior;
     a=uq_audit_outputs(out); flags=flags+a.flags; initialflags=initialflags+a.initial_flags; fallbacks=fallbacks+a.fallbacks; maxres=max(maxres,a.max_normal);
     for route=1:2
      if route==1, S=out.direct; tag='Direct'; else, S=out.split; tag='Split'; end
      Z(:,:,route,ip)=S; row=uq_metric_row(S,scene,c,c.model_name,suite,seed,methods{ip},tag,out,prior,preset,est);
      row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma); count=count+1;
      if isempty(rows), rows=row; else, rows(count)=row; end
     end
     ev=struct('Model',c.model_name,'Suite',suite,'Scenario',name,'Seed',seed,'Method',methods{ip},'Preset',preset.name, ...
      'ScaleGrid',out.scale_grid,'ResidualGrid',out.residual_grid,'Target',out.target,'SelectedScale',out.scale,'Audit',a);
     if isempty(evidence), evidence=ev; else, evidence(end+1)=ev; end
     maxfit=max(maxfit,out.residual_mismatch); endpoints=endpoints+out.endpoint;
    end
    if isu==1 && seed==o.test_seeds(1)
     uq_writejson(fullfile(outdir,['example_' c.model_name '_' name '.json']),struct('Model',c.model_name,'Scenario',name, ...
      't',c.t,'offsets',scene.offsets,'truth',scene.truth,'observed',D,'Z',Z,'methods',{methods},'estimated_Q',est.effective_Q, ...
      'estimated_frequency',est.estimated_frequency,'public_config',p,'presets',presets, ...
      'lambda',prior.lambda,'estimated_sigma',prior.sigma));
    end
    inputs=inputs+1; writetable(struct2table(rows),fullfile(outdir,'metrics_partial.csv'));
    uq_writejson(fullfile(outdir,'run_status.json'),struct('status','RUNNING','rows',count,'inputs',inputs,'model',c.model_name,'suite',suite,'scenario',name,'seed',seed,'elapsed_seconds',toc(start)));
    fprintf('Test %s %s %s seed %d: %d rows %.1fs\n',c.model_name,suite,name,seed,count,toc(start));
   end
  end
 end
end
writetable(struct2table(rows),fullfile(outdir,'metrics.csv')); uq_writejson(fullfile(outdir,'selection_evidence.json'),evidence);
status=struct('status','COMPLETED','rows',count,'inputs',inputs,'solver_nonzero_flags',flags,'initial_nonzero_flags',initialflags, ...
 'fallback_count',fallbacks,'max_normal_residual',maxres,'max_target_mismatch',maxfit,'scale_endpoints',endpoints,'elapsed_seconds',toc(start));
uq_writejson(fullfile(outdir,'run_status.json'),status); assert(count==o.expected_test_rows && flags==0 && maxres<2e-7);
end
