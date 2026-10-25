function p1_run_grid(outdir,spec,smoke)
mkdir(outdir); start=tic;
uq_writejson(fullfile(outdir,'executed_spec.json'),spec);
splits={'historical','development','heldout_new_geology'};
rows=struct([]); count=0; flags=0; maxnormal=0; audits=struct([]); na=0;
for is=1:numel(splits)
 sname=splits{is}; listed=spec.models.(sname);
 if isempty(listed), continue; end
 [configs,~]=p1_models(sname);
 for im=1:numel(listed)
  c=configs{im}; public=uq_public_config(c,c.nominal_frequency);
  trueModels=prepare_models_v2(c.offset,c);
  [~,~,~,calScene]=uq_case('SmoothAVO',c,trueModels);
  previousModels=[];
  for ic=1:numel(spec.calibration_seeds)
   cseed=spec.calibration_seeds(ic);
   [calD,~]=uq_observation(calScene,c,trueModels,1:11,cseed,18);
   [est,cm,invmodels,R_avo]=p1_prepare_calibration(calD,c.offset,public);
   record=struct('model',c.model_name,'calibration_seed',cseed,'Q',est.effective_Q, ...
    'frequency',est.estimated_frequency,'operator_change_from_previous',0);
   if ~isempty(previousModels)
    record.operator_change_from_previous=norm(invmodels{6}.G-previousModels{6}.G,'fro');
   end
   previousModels=invmodels;
   uq_writejson(fullfile(outdir,sprintf('calibration_%s_%d.json',c.model_name,cseed)),est);
   na=na+1; if na==1, audits=record; else, audits(na)=record; end
   ns=2+numel(spec.stress_suites); if smoke, ns=1+numel(spec.stress_suites); end
   for isu=1:ns
    suite='Main'; snr=18; scenarios=spec.scenarios; stress=[];
    cm_use=cm; models_use=invmodels; R_use=R_avo; est_use=est;
    if ~smoke && isu==2
     suite='LowerSNR'; snr=12; scenarios=spec.stress_scenarios;
    elseif (smoke && isu>1) || (~smoke && isu>2)
     idx=isu-2; if smoke, idx=isu-1; end
     stress=spec.stress_suites(idx); suite=stress.name;
     scenarios=spec.stress_scenarios; if smoke, scenarios=scenarios(1); end
     if startsWith(suite,'StressVelocity')
      bad=public; bad.vp=public.vp*stress.vp_multiplier;
      [est_use,cm_use,models_use,R_use]=p1_prepare_calibration(calD,c.offset,bad);
      assert(norm(cm_use.vp-cm.vp)>0 && norm(models_use{6}.G-invmodels{6}.G,'fro')>1e-8);
      uq_writejson(fullfile(outdir,sprintf('velocity_stress_%s_%s_%d.json',stress.name,c.model_name,cseed)), ...
       struct('vp_multiplier',stress.vp_multiplier,'baseline_vp',cm.vp,'stress_vp',cm_use.vp, ...
       'operator_frobenius_change',norm(models_use{6}.G-invmodels{6}.G,'fro'),'calibration',est_use));
     end
    end
    if ~strcmp(sname,'historical'), suite=[sname '_' suite]; end
    for isc=1:numel(scenarios)
     scname=char(scenarios{isc});
     [models,domain,map,scene]=uq_case(scname,c,models_use);
     for seed=reshape(spec.test_seeds,1,[])
      [D,trueSigma]=uq_observation(scene,c,trueModels,map,seed,snr,stress);
      for ib=1:numel(spec.methods)
       bname=spec.methods(ib).name;
       [out,prior,preset]=p1_solve_branch(D,models,map,cm_use,domain,R_use,bname);
       au=uq_audit_outputs(out); flags=flags+au.flags; maxnormal=max(maxnormal,au.max_normal);
       for ir=1:2
        route='Direct'; S=out.direct; if ir==2, route='Split'; S=out.split; end
        row=uq_metric_row(S,scene,c,c.model_name,suite,seed,bname,route,out,prior,preset,est_use);
        row.EstimatedSigmaRatio=mean(prior.sigma./trueSigma); row.CalibrationSeed=cseed;
        count=count+1; if count==1, rows=row; else, rows(count)=row; end
       end
       fprintf('%s %s %s seed=%d cal=%d error=%.9f normal=%.3e/%.3e elapsed=%.1fs\n', ...
        c.model_name,suite,[scname '/' bname],seed,cseed,row.RelativeError, ...
        out.direct_record.normal_residual,out.split_record.normal_residual,toc(start));
      end
      writetable(struct2table(rows),fullfile(outdir,'metrics.csv'));
     end
    end
   end
  end
 end
end
uq_writejson(fullfile(outdir,'calibration_operator_audit.json'),audits);
uq_writejson(fullfile(outdir,'checks.json'),struct('solver_flags',flags,'max_normal_residual',maxnormal, ...
 'rows',count,'elapsed_seconds',toc(start),'smoke',smoke,'version',version));
assert(flags==0 && maxnormal<=2e-7,'Solver quality failure; outputs retained for diagnosis');
fprintf('COMPLETED %d rows %.1fs\n',count,toc(start));
end
