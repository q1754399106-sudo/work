function probe_baselines()
% Smoke test of SBL-TVD and Wang stabilized inverse Q on one development block.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(2);
soft=fullfile(fileparts(root),'soft_refit_20261023');
for which={'primary','complex'}
 [cm,ref]=load_dev(soft,which{1});data=ref.datasets(1);ids=1:11;D=data.D(:,ids);T=data.target(:,ids);sig=reshape(data.sigma(ids),1,[]);
 used=cm;used.Q=reshape(ref.reference_Q,1,[]);d=prune_dictionary(used,450);
 tic;full=reflectivity_refit_prepared(D,repmat({d},1,11),sig);t1=toc;
 tic;sbl=sbl_tvd(D,d,sig);t2=toc;
 e=@(A)norm(A-T,'fro')/norm(T,'fro');
 fprintf('%s Full %.3f (%.1fs) SBL %.3f (%.1fs, active %d, iters max %d)\n',which{1},e(full.gather),t1,e(sbl.gather),t2,max(sbl.active),max(sbl.iterations));
 for G=[10 20 30 40],fprintf('  Wang G=%d: %.3f\n',G,e(wang_inverse_q(D,used,450,G)));end
 fprintf('  raw data vs target: %.3f\n',e(D));
end
end
