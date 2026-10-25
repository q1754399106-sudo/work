function run_complex_model(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
input=jsondecode(fileread(fullfile(root,'inputs','complex_input.json')));p=input.public;
for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end
p.t=p.t(:);assert(~any(isfield(p,{'Q','fpeak','truth','rho','vs'})));
assert(~isfolder(outdir));mkdir(outdir);start=tic;
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));
[est,corrected,cal]=p1_prepare_candidate(input.calibrationD,p.offset,p,policy);
cmfull=uq_inverse_config(corrected,est);
estbase=p1_candidate_domain_estimate(cal.fixed_estimate,est,corrected,policy,'Section');cmbase=uq_inverse_config(corrected,estbase);
offsets=repmat(input.section_offset_m,1,11);
op=fixed_v1_operators(cmbase,offsets);shared=fixed_v1_operators(cmfull,offsets);
dict=prune_dictionary(cmfull,input.section_offset_m);dicts=repmat({dict},1,11);
uq_writejson(fullfile(outdir,'calibration.json'),struct('estimate',est,'original_estimate',estbase,'corrected',corrected,'diagnostics',cal));
fprintf('Calibration/dictionaries ready elapsed=%.1fs\n',toc(start));
total=0;flags=0;maxnormal=0;
for id=1:numel(input.datasets)
 data=input.datasets(id);D=data.D;[nt,nx]=size(D);assert(nt==p.nt&&mod(nx,11)==0);
 fixed=zeros(size(D));fixed_shared=fixed;unpruned=fixed;optimized=fixed;diagnostics=cell(1,nx/11);
 for k=1:nx/11
  idx=(k-1)*11+(1:11);block=D(:,idx);sigma=zeros(1,11);
  for j=1:11,ni=estimate_trace_snr(block(:,j),p.t);sigma(j)=ni.sigma;end
  [fixed(:,idx),~,~]=fixed_v1_apply(block,op);[fixed_shared(:,idx),~,~]=fixed_v1_apply(block,shared);
  initial=reflectivity_refit_prepared(block,dicts,sigma);unpruned(:,idx)=initial.gather;
  out=reflectivity_prune_prepared(block,dicts,sigma,initial);optimized(:,idx)=out.gather;
  flags=flags+out.optimization_flags;maxnormal=max(maxnormal,out.normal_residual);
  diagnostics{k}=struct('block',k,'traces',idx,'sigma',sigma,'flags',out.optimization_flags,'normal',out.normal_residual,'support_before',out.support_before_pruning,'support_after',out.support_after,'removed',out.removed_coefficients);
  total=total+1;uq_writejson(fullfile(outdir,'progress.json'),struct('blocks',total,'expected',numel(input.datasets)*nx/11,'dataset',data.name,'elapsed_seconds',toc(start)));
  fprintf('%s block %d/%d total=%d flags=%d elapsed=%.1fs\n',data.name,k,nx/11,total,out.optimization_flags,toc(start));
 end
 observed=D;t=p.t;
 save(fullfile(outdir,[data.name '.mat']),'observed','fixed','fixed_shared','unpruned','optimized','t','diagnostics','-v7');
 uq_writejson(fullfile(outdir,[data.name '.json']),struct('name',data.name,'t',t,'observed',observed,'fixed',fixed,'fixed_shared',fixed_shared,'unpruned',unpruned,'optimized',optimized,'diagnostics',{diagnostics}));
end
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'blocks',total,'optimization_flags',flags,'max_normal',maxnormal,'elapsed_seconds',toc(start)));
end
