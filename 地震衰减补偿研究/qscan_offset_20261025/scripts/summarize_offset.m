function summarize_offset()
% Scoring (PLAN_ZH.md): per model and scene, the Standard error curve over the 14 factors (R=6),
% QScan + Standard (D), D moved down by k=1 and k=2 grid steps (S1, S2), nominal Standard (B),
% the grid offset between QScan's modal choice and the error-optimal grid factor, the noise-free
% optimum and the strong-reflection gain; judgments O1-O3.
P=os_setup();A=fullfile(P.root,'analysis');if ~isfolder(A),mkdir(A);end
f=P.factors;nf=numel(f);R=P.R;step=log(P.grid(2)/P.grid(1));curve={};summ={};
for m=1:numel(P.models)
 Tt=load(fullfile(P.root,'truth',[P.models(m).name '_truth.mat']),'target');T=Tt.target;ET=sum(T(:).^2);
 Lt=sqrt(movmean(T.^2,21,1));strong=Lt>=.25*max(Lt(:));ES=sum(T(strong).^2);
 for s=1:numel(P.scenes)
  sc=P.scenes(s);name=[P.models(m).name '_' sc.tag];dir_=fullfile(P.root,'cache','sweep',name);
  sq=zeros(R,nf);sqA=sq;gain=sq;best=zeros(R,1);
  for r=1:R
   L=load(fullfile(dir_,sprintf('r%03d.mat',r)));assert(L.info.optimization_flags==0);best(r)=L.info.best;
   for k=1:nf
    X=double(L.STD(:,:,k));sq(r,k)=sum((X-T).^2,'all')/ET;gain(r,k)=sum(X.*T.*strong,'all')/ES;
    sqA(r,k)=sum((double(L.FCR(:,:,k))-T).^2,'all')/ET;
   end
  end
  N=load(fullfile(dir_,'r000.mat'));nfree=zeros(1,nf);for k=1:nf,nfree(k)=sum((double(N.STD(:,:,k))-T).^2,'all')/ET;end
  mc=mean(sq,1);q=f*sc.q_scale;
  for k=1:nf
   curve(end+1,:)={P.models(m).name,sc.tag,k,f(k),q(k),mc(k),std(sq(:,k))/sqrt(R),mean(sqA(:,k)),nfree(k),mean(gain(:,k)),sum(best==k-1)}; %#ok<AGROW>
  end
  rows=@(idx) sub2ind([R nf],(1:R)',1+idx);
  mD=mean(sq(rows(best)));m1=mean(sq(rows(max(best-1,1))));m2=mean(sq(rows(max(best-2,1))));mB=mc(1);
  [mopt,io]=min(mc);[~,ig]=min(mc(2:end));qscan=mode(best);[~,inf_]=min(nfree);
  summ(end+1,:)={P.models(m).name,sc.tag,mD,m1,m2,mB,mopt,m1/mD-1,m2/mD-1,mopt/mD-1,mB/mD-1, ...
   q(io),P.grid(ig)*sc.q_scale,P.grid(qscan)*sc.q_scale,qscan-ig,sum(best==qscan),q(inf_), ...
   mean(gain(rows(best))),mean(gain(rows(max(best-1,1)))),mean(gain(rows(max(best-2,1)))),mean(gain(:,io))}; %#ok<AGROW>
  fprintf('%s: D %.4f S1 %.4f S2 %.4f B %.4f opt %.4f | q* %.3f q_scan %.3f offset %d steps | noise-free q* %.3f\n', ...
   name,mD,m1,m2,mB,mopt,q(io),P.grid(qscan)*sc.q_scale,qscan-ig,q(inf_));
 end
end
writetable(cell2table(curve,'VariableNames',{'model','scene','index','factor','q_ratio_to_Q0','mse','mse_se','mse_fullcalrefit', ...
 'noise_free_mse','strong_gain','times_selected'}),fullfile(A,'offset_curves.csv'));
S=cell2table(summ,'VariableNames',{'model','scene','mse_D','mse_S1','mse_S2','mse_B','mse_opt','rel_S1','rel_S2','rel_opt','rel_B', ...
 'q_best','q_best_grid','q_qscan_mode','offset_steps','qscan_mode_count','q_best_noise_free','gain_D','gain_S1','gain_S2','gain_opt'});
writetable(S,fullfile(A,'offset_summary.csv'));
p=S(strcmp(S.scene,'Q100_14dB'),:);k=p.offset_steps;out=struct();
[cnt,val]=groupcounts(k);[top,i]=max(cnt);
out.O1=struct('offset_steps',k','modal_offset',val(i),'modal_count',top,'spread',max(k)-min(k), ...
 'verdict',ternary(max(k)-min(k)<=1&&top>=6,'offset stable across models','offset not stable across models'));
for name={'S1','S2'}
 rel=S.(['rel_' name{1}]);ok=sum(rel<0)>=20&&max(rel)<=.05;
 out.O2.(name{1})=struct('n_improved',sum(rel<0),'n',numel(rel),'worst_rel',max(rel),'median_rel',median(rel), ...
  'verdict',ternary(ok,'safe shift','not safe'));
end
near=abs(log(p.q_best_noise_free))<=step;
out.O3=struct('q_best_noise_free',p.q_best_noise_free','n_near_Q0',sum(near),'verdict', ...
 ternary(sum(near)>=6,'noise-free optimum at the upscaled Q0','noise-free optimum not at the upscaled Q0'));
fid=fopen(fullfile(A,'decisions.json'),'w');fprintf(fid,'%s',jsonencode(out,'PrettyPrint',true));fclose(fid);disp(jsonencode(out,'PrettyPrint',true));
end

function v=ternary(c,a,b)
if c,v=a;else,v=b;end
end
