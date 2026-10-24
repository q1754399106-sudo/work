function probe_feasibility()
% Development probe on the archived F1/F2 Q0 observations (already-seen data): does the
% Q0-world bootstrap give the sign of the replicate-truth Delta MSE? Also checks, with truth,
% the premise that FullCalRefit's bias is nearly orthogonal to the QScan change at Q0.
P=bd_setup();maxNumCompThreads(P.threads);rows={};
dec=readtable(fullfile(P.in,'qscan_bias_variance_20261024','analysis','decomposition.csv'),'TextType','string');
mom=fullfile(fileparts(P.root),'qscan_bias_variance_20261024','cache');
for model=["F1","F2"]
 for cond=["Moderate","Strong"]
  sc=load_obs(P,'bae',char(model),char(cond));t=sc.cm.t;tic;
  x1=full_refit(sc.D,sc.sigma,sc.cm,sc.Q0,sc.offset);sel=qscan_select(sc.D,sc.sigma,sc.cm,sc.Q0,sc.offset);
  c=P.grid(sel.best);xs=full_refit(sc.D,sc.sigma,sc.cm,sc.Q0*c,sc.offset);phi=noise_ar1(sc.D,t,P.window);
  seeds=20261025e2+(1:P.K);
  H0=bootstrap_world(x1,sc.Q0,sc.Q0,sc.Q0*c,sc.cm,sc.offset,phi,sc.sigma,seeds);
  H1=bootstrap_world(xs,sc.Q0*c,sc.Q0,sc.Q0*c,sc.cm,sc.offset,phi,sc.sigma,seeds);
  d=dec(dec.model==model&dec.condition==cond,:);mb=d.mse_rel(d.estimator=="c=1.0000");ms=d.mse_rel(d.estimator=="QScan");
  b2=[d.bias2_rel(d.estimator=="c=1.0000"),d.bias2_rel(d.estimator=="QScan")];
  T=load(fullfile(P.in,'bae_q_20261024','truth','independent',[char(model) '_truth.mat']));T=T.tr(strcmp({T.tr.name},cond)).target;
  M=load(fullfile(mom,sprintf('moments_%s_%s.mat',model,cond)),'Xbar');Xb=double(M.Xbar);dm=Xb(:,:,end)-Xb(:,:,1);b1=Xb(:,:,1)-T;ET=sum(T(:).^2);
  rows(end+1,:)={char(model),char(cond),c,phi,(ms-mb)/mb,b2(2)-b2(1),sum(dm(:).^2)/ET,2*sum(b1(:).*dm(:))/ET, ...
   H0.rel_delta_mse,H0.bias2(2)-H0.bias2(1),H0.var(2)-H0.var(1),H1.rel_delta_mse,toc}; %#ok<AGROW>
  fprintf('%s %s c=%.4f phi=%.3f | truth dMSE %+.3f (dBias2 %+.4f = |dm|^2 %.4f + 2b.dm %+.4f) | H0 %+.3f (dB2 %+.4f dVar %+.4f) | H1 %+.3f | %.0fs\n',rows{end,:});
 end
end
writetable(cell2table(rows,'VariableNames',{'model','condition','qscan_factor','phi_hat','truth_rel_delta_mse','truth_delta_bias2', ...
 'truth_mean_change_sq','truth_cross_term','H0_rel_delta_mse','H0_delta_bias2','H0_delta_var','H1_rel_delta_mse','seconds'}), ...
 fullfile(P.root,'analysis','probe_feasibility.csv'));
fprintf('DONE probe\n');
end
