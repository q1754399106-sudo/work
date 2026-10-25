function summarize_fd()
% Scoring (PLAN_ZH.md): per model and scene, the Standard and FullCalRefit error curves over the
% 14 factors (R=6), QScan + Standard (D), D moved down by 1 and 2 grid steps (S1, S2), nominal
% Standard (B), nominal FullCalRefit (A), the grid offset between QScan's modal choice and the
% error-optimal grid factor, the noise-free optimum and the strong-reflection gain; judgments
% E1-E4. Statement for statement as qscan_offset_20261025/scripts/summarize_offset.m where the
% quantities coincide (fd_setup, six models, A and E3/E4 added).
P=fd_setup();A=fullfile(P.root,'analysis');if ~isfolder(A),mkdir(A);end
f=P.factors;nf=numel(f);R=P.R;step=log(P.grid(2)/P.grid(1));curve={};summ={};nm=numel(P.models);
for m=1:nm
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
  mc=mean(sq,1);mA=mean(sqA,1);q=f*sc.q_scale;
  for k=1:nf
   curve(end+1,:)={P.models(m).name,sc.tag,k,f(k),q(k),mc(k),std(sq(:,k))/sqrt(R),mA(k),nfree(k),mean(gain(:,k)),sum(best==k-1)}; %#ok<AGROW>
  end
  rows=@(idx) sub2ind([R nf],(1:R)',1+idx);
  mD=mean(sq(rows(best)));m1=mean(sq(rows(max(best-1,1))));m2=mean(sq(rows(max(best-2,1))));mB=mc(1);
  [mopt,io]=min(mc);[~,ig]=min(mc(2:end));qscan=mode(best);[~,inf_]=min(nfree);[~,ioA]=min(mA);
  summ(end+1,:)={P.models(m).name,sc.tag,mD,m1,m2,mB,mA(1),mopt,m1/mD-1,m2/mD-1,mopt/mD-1,mB/mD-1,mB/mA(1)-1, ...
   q(io),P.grid(ig)*sc.q_scale,P.grid(qscan)*sc.q_scale,qscan-ig,sum(best==qscan),q(inf_),q(ioA), ...
   mean(gain(rows(best))),mean(gain(rows(max(best-1,1)))),mean(gain(rows(max(best-2,1)))),mean(gain(:,io))}; %#ok<AGROW>
  fprintf('%s: D %.4f S1 %.4f S2 %.4f B %.4f A %.4f opt %.4f | q* %.3f q_scan %.3f offset %d steps | noise-free q* %.3f | FullCalRefit q* %.3f\n', ...
   name,mD,m1,m2,mB,mA(1),mopt,q(io),P.grid(qscan)*sc.q_scale,qscan-ig,q(inf_),q(ioA));
 end
end
writetable(cell2table(curve,'VariableNames',{'model','scene','index','factor','q_ratio_to_Q0','mse','mse_se','mse_fullcalrefit', ...
 'noise_free_mse','strong_gain','times_selected'}),fullfile(A,'fd_curves.csv'));
S=cell2table(summ,'VariableNames',{'model','scene','mse_D','mse_S1','mse_S2','mse_B','mse_A','mse_opt','rel_S1','rel_S2','rel_opt','rel_B','rel_B_vs_A', ...
 'q_best','q_best_grid','q_qscan_mode','offset_steps','qscan_mode_count','q_best_noise_free','q_best_fullcalrefit','gain_D','gain_S1','gain_S2','gain_opt'});
writetable(S,fullfile(A,'fd_summary.csv'));
p=S(strcmp(S.scene,'Q100_14dB'),:);u=S(strcmp(S.scene,'Q070_14dB'),:);k=p.offset_steps;out=struct();
[cnt,val]=groupcounts(k);[top,i]=max(cnt);
out.E1=struct('offset_steps',k','modal_offset',val(i),'modal_count',top,'spread',max(k)-min(k),'both_signs',any(k>0)&&any(k<0), ...
 'verdict',ternary(max(k)-min(k)<=1&&top>=5,'offset stable across models','offset not stable across models'));
near=abs(log(p.q_best_noise_free))<=step;
out.E2=struct('q_best_noise_free',p.q_best_noise_free','n_near_Q0',sum(near),'verdict', ...
 ternary(sum(near)>=5,'noise-free optimum at the upscaled Q0','noise-free optimum not at the upscaled Q0'));
better=p.mse_B<p.mse_A;
out.E3=struct('rel_B_vs_A',p.rel_B_vs_A','n_standard_better',sum(better),'verdict', ...
 ternary(sum(better)>=5,'coherence shrinkage beats FullCalRefit at nominal Q','coherence shrinkage does not beat FullCalRefit at nominal Q'));
bQ0=p.mse_B<p.mse_D;dQ7=u.mse_D<u.mse_B;
out.E4=struct('n_B_better_Q0',sum(bQ0),'rel_B_Q0',p.rel_B','n_D_better_Q070',sum(dQ7),'rel_B_Q070',u.rel_B', ...
 'verdict_Q0',ternary(sum(bQ0)>=4,'nominal Standard better when Q is correct','nominal Standard not better when Q is correct'), ...
 'verdict_Q070',ternary(sum(dQ7)>=5,'QScan + Standard better when Q is 30% low','QScan + Standard not better when Q is 30% low'));
for name={'S1','S2'}
 rel=S.(['rel_' name{1}]);
 out.shifts.(name{1})=struct('n_improved',sum(rel<0),'n',numel(rel),'worst_rel',max(rel),'median_rel',median(rel));
end
fid=fopen(fullfile(A,'decisions.json'),'w');fprintf(fid,'%s',jsonencode(out,'PrettyPrint',true));fclose(fid);disp(jsonencode(out,'PrettyPrint',true));
end

function v=ternary(c,a,b)
if c,v=a;else,v=b;end
end
