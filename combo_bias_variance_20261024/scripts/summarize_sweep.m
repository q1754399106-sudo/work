function summarize_sweep()
% Stage 2 scoring (PLAN_STAGE2_ZH.md): per scene and factor, MSE = bias^2 + variance of the
% Standard and FullCalRefit outputs (R=20) and strong-reflection gain; QScan + Standard (D);
% the cross-fitted best fixed factor; in-sample and per-replicate oracles; the noise-free curve;
% K1-K3.
P=cb_setup();A=fullfile(P.root,'analysis');factors=[1 P.grid];nf=numel(factors);R=20;half={1:10,11:20};
curve={};summ={};
for s=[1 2 5 6 3 4 7 8]
 sc=load_scene(P,s);T=sc.T;ET=sum(T(:).^2);Lt=sqrt(movmean(T.^2,21,1));strong=Lt>=.25*max(Lt(:));ES=sum(T(strong).^2);
 [n,nx]=size(T);X=zeros(n,nx,nf,R,'single');Y=X;best=zeros(R,1);J=zeros(R,numel(P.grid));
 for r=1:R
  L=load(fullfile(P.root,'cache','sweep',sc.name,sprintf('r%03d.mat',r)));assert(L.info.stage1_identical==1);
  X(:,:,:,r)=L.STD;Y(:,:,:,r)=L.FCR;best(r)=L.info.best;J(r,:)=L.info.objective;
 end
 N=load(fullfile(P.root,'cache','sweep',sc.name,'r000.mat'));
 sq=zeros(R,nf);sqA=sq;
 for f=1:nf
  for k=1:2
   if k==1,Z=double(squeeze(X(:,:,f,:)));lab='Standard';else,Z=double(squeeze(Y(:,:,f,:)));lab='FullCalRefit';end
   e=squeeze(sum((Z-T).^2,[1 2]))/ET;if k==1,sq(:,f)=e;else,sqA(:,f)=e;end
   Zb=mean(Z,3);va=sum(var(Z,0,3),'all')/ET;b2=sum((Zb-T).^2,'all')/ET-va/R;
   gm=sum(Zb.*T.*strong,'all')/ES;nfree=sum((double(N.STD(:,:,f))-T).^2,'all')/ET;if k==2,nfree=sum((double(N.FCR(:,:,f))-T).^2,'all')/ET;end
   curve(end+1,:)={sc.name,lab,f,factors(f),factors(f)*sc.q_scale,mean(e),b2,va,gm,nfree,f>1&&any(best==f-1),sum(best==f-1)}; %#ok<AGROW>
  end
 end
 rows=sub2ind([R nf],(1:R)',1+best);mD=mean(sq(rows));
 [~,fin]=min(mean(sq,1));
 cf=zeros(1,2);eval_=[];for k=1:2,[~,cf(k)]=min(mean(sq(half{k},:),1));eval_=[eval_;sq(half{3-k},cf(k))];end %#ok<AGROW>
 mfix=mean(eval_);morc=mean(min(sq,[],2));[~,fnf]=min(cellfun(@(c)c,curve(strcmp(curve(:,1),sc.name)&strcmp(curve(:,2),'Standard'),10)));
 [~,fA]=min(mean(sqA,1));[~,jbest]=min(mean(J,1));qscan=mode(best);
 summ(end+1,:)={sc.name,ismember(s,[1 2 5 6]),mD,mfix,1-mfix/mD,mean(sq(:,fin)),1-mean(sq(:,fin))/mD,morc,1-morc/mD, ...
  factors(cf(1))*sc.q_scale,factors(cf(2))*sc.q_scale,factors(fin)*sc.q_scale,P.grid(qscan)*sc.q_scale,factors(fnf)*sc.q_scale, ...
  factors(fA)*sc.q_scale,P.grid(jbest)*sc.q_scale}; %#ok<AGROW>
 fprintf('%s: D %.4f fixed(cf) %.4f G %.3f | q* %.3f/%.3f in-sample %.3f qscan %.3f noise-free %.3f FullCal %.3f\n',sc.name,mD,mfix,1-mfix/mD, ...
  factors(cf(1))*sc.q_scale,factors(cf(2))*sc.q_scale,factors(fin)*sc.q_scale,P.grid(qscan)*sc.q_scale,factors(fnf)*sc.q_scale,factors(fA)*sc.q_scale);
end
writetable(cell2table(curve,'VariableNames',{'scene','method','index','factor','q_ratio_to_Q0','mse','bias2','var','strong_gain', ...
 'noise_free_mse','selected_by_qscan','times_selected'}),fullfile(A,'sweep_curves.csv'));
S=cell2table(summ,'VariableNames',{'scene','primary','mse_D','mse_fixed_crossfit','G_crossfit','mse_fixed_insample','G_insample', ...
 'mse_oracle_per_replicate','G_oracle','q_best_half1','q_best_half2','q_best_insample','q_qscan_mode','q_best_noise_free', ...
 'q_best_fullcalrefit','q_min_objective'});
writetable(S,fullfile(A,'sweep_summary.csv'));
p=S(S.primary,:);out=struct();
out.K1=struct('G_crossfit',p.G_crossfit','mean_G',mean(p.G_crossfit),'verdict',ternary(mean(p.G_crossfit)>=.05,'selection mechanism worth studying','limited room; stop'));
lower=(p.q_best_half1<p.q_qscan_mode)+(p.q_best_half2<p.q_qscan_mode);
out.K2=struct('q_best_halves',[p.q_best_half1 p.q_best_half2],'q_qscan',p.q_qscan_mode','verdict', ...
 ternary(sum(lower==2)>=3,'QScan systematically selects too high a Q (under-compensation)','not established'));
fid=fopen(fullfile(A,'stage2_decisions.json'),'w');fprintf(fid,'%s',jsonencode(out,'PrettyPrint',true));fclose(fid);disp(jsonencode(out,'PrettyPrint',true));
end

function v=ternary(c,a,b)
if c,v=a;else,v=b;end
end
