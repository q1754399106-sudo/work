function time_probe()
% Timing of one 11-trace block per model (development data), no scoring.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
test_bae_math();soft=fullfile(fileparts(root),'soft_refit_20261023');
for which={'primary','complex'}
 [cm,ref]=load_dev(soft,which{1});data=ref.datasets(2);ids=1:11;D=data.D(:,ids);sig=reshape(data.sigma(ids),1,[]);
 used=cm;used.Q=reshape(ref.reference_Q,1,[])*0.7;d=prune_dictionary(used,450);
 tic;initial=reflectivity_refit_prepared(D,repmat({d},1,11),sig);t1=toc;
 tic;Gs=q_sample_operators(used,450,.2,32,.7,20261024);t2=toc;
 tic;b=bae_q_refit(D,d,Gs,sig,initial,true);t3=toc;
 fprintf('%s nt=%d refit %.1fs operators %.1fs bae %.1fs rank %d flags %d\n',which{1},cm.nt,t1,t2,t3,max(b.error_rank),b.optimization_flags);
end
end
