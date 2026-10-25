function timing_test()
% Development timing of the diagnostic parts on the archived F1/Moderate gather (machine shared
% with other jobs, so absolute times are inflated; the ratios matter).
P=bd_setup();maxNumCompThreads(P.threads);sc=load_obs(P,'bae','F1','Moderate');cm=sc.cm;
tic;sel=qscan_select(sc.D,sc.sigma,cm,sc.Q0,sc.offset);t_sel=toc;
tic;x1=full_refit(sc.D,sc.sigma,cm,sc.Q0,sc.offset);t_ref=toc;
fold=mod(ceil((1:cm.nt)'/5)-1,5)+1;obs=fold~=1;
c=cm;c.Q=sc.Q0;d=prune_dictionary(c,sc.offset);B=d.B(obs,:);Gram=B'*B;Gram=(Gram+Gram')/2;L=norm(Gram,1);
lambda=.5*sqrt(2*log(cm.nt))*sc.sigma;tic;[~,di]=refit_lasso(Gram,B'*sc.D(obs,1:11),lambda(1:11),L);t_lasso=toc;
tic;masked_refit(sc.D,sc.sigma,cm,sc.Q0,sc.offset,obs);t_mask=toc;
fprintf('select %.1fs, refit %.1fs, masked lasso block %.1fs (iterations %s, fallback %d), masked refit %.1fs\n', ...
 t_sel,t_ref,t_lasso,mat2str(di.iterations),sum(di.coordinate_fallback),t_mask);
end
