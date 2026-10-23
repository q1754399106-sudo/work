function run_soft_trial()
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));maxNumCompThreads(4);
test_soft_math();
old=fullfile(fileparts(root),'q_dense_gathers_20261023');started=tic;total=0;flags=0;maxdiff=0;maxnormal=0;
for protocol={'primary','complex'}
 which=protocol{1};input=jsondecode(fileread(fullfile(root,'inputs',[which '.json'])));
 ref=jsondecode(fileread(fullfile(old,'results',which,'reference.json')));datasets=ref.datasets;
 if strcmp(which,'primary')
  cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';
 else
  p=input.public;for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end
  p.t=p.t(:);cm=uq_inverse_config(p,input.estimate);
 end
 referenceQ=reshape(ref.reference_Q,1,[]);outdir=fullfile(root,'results',which);
 if ~isfolder(outdir),mkdir(outdir);end
 uq_writejson(fullfile(outdir,'reference.json'),ref);
 for errorpct=reshape(input.errors_percent,1,[])
  used=cm;used.Q=referenceQ*(1+errorpct/100);dict=prune_dictionary(used,450);
  Wscaled=dict.W./reshape(dict.scale,1,[]);dict.source_gram=Wscaled'*Wscaled;dicts=repmat({dict},1,11);
  baseline=jsondecode(fileread(fullfile(old,'results',which,sprintf('Q_%+d.json',errorpct))));
  records=cell(1,numel(datasets));
  for k=1:numel(datasets)
   data=datasets(k);D=data.D;sig=reshape(data.sigma,1,[]);assert(size(D,2)==99);
   prior=baseline.records(k);assert(strcmp(prior.name,data.name));
   soft=zeros(size(D));unpruned=soft;diagnostics=cell(1,9);
   for block=1:9
    ids=(block-1)*11+(1:11);cache=fullfile(root,'cache',sprintf('%s_Q%+d_%s_%02d.mat',which,errorpct,data.name,block));
    if isfile(cache),cached=load(cache,'initial','out');initial=cached.initial;out=cached.out;
    else
     initial=reflectivity_refit_prepared(D(:,ids),dicts,sig(ids));
     out=reflectivity_soft_prepared(D(:,ids),dicts,sig(ids),initial);
     save(cache,'initial','out','-v7');
    end
    difference=norm(initial.gather-prior.unpruned(:,ids),'fro')/max(norm(prior.unpruned(:,ids),'fro'),realmin);
    assert(difference<1e-6,'Recomputed FullCalRefit differs from frozen baseline.');
    soft(:,ids)=out.gather;unpruned(:,ids)=initial.gather;
    flags=flags+initial.optimization_flags;total=total+1;maxdiff=max(maxdiff,difference);maxnormal=max(maxnormal,out.initial_normal_residual);
    diagnostics{block}=struct('group_gain',{out.group_gain},'group_energy',{out.group_energy},'group_noise_energy',{out.group_noise_energy}, ...
       'cached_baseline_relative_difference',difference,'initial_normal',out.initial_normal_residual);
    progress=struct('protocol',which,'error_percent',errorpct,'condition',data.name,'block',block,'blocks_done',total,'blocks_expected',180, ...
       'flags',flags,'max_baseline_difference',maxdiff,'elapsed_seconds',toc(started));
    uq_writejson(fullfile(root,'progress.json'),progress);
    fprintf('%s Q%+d %s %d/9 total=%d/180 flags=%d diff=%.3g elapsed=%.1fs\n',which,errorpct,data.name,block,total,flags,maxdiff,toc(started));
   end
   records{k}=struct('name',data.name,'fixed',prior.fixed,'unpruned',unpruned,'hard',prior.optimized,'soft',soft,'diagnostics',{diagnostics});
  end
  uq_writejson(fullfile(outdir,sprintf('Q_%+d.json',errorpct)),struct('error_percent',errorpct,'used_Q',used.Q,'records',{records}));
 end
end
checks=struct('complete',true,'blocks',total,'expected_blocks',180,'optimization_flags',flags,'max_cached_baseline_difference',maxdiff,'max_initial_normal',maxnormal,'elapsed_seconds',toc(started));
uq_writejson(fullfile(root,'checks.json'),checks);
uq_writejson(fullfile(root,'RUN_LEDGER.json'),struct('status','computed','completed_blocks',total,'expected_blocks',180,'outputs_root',root));
end
