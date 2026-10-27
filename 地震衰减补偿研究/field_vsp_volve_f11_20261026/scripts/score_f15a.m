function score_f15a()
% F-15A addition (PLAN_ZH.md "F-15A 补算", descriptive): T1/T2 of SC and QScan+SC on the experiment-5 curtain,
% scored exactly as experiment 5 (its archived references.mat, curtain idx and 6-70 Hz filter). The archived
% Input/Standard/QScan_Standard outputs are rescored first and must reproduce experiment 5's scores.csv.
P=setup_paths();root=P.root;E=P.f15a;in=load(fullfile(E,'inputs','curtain.mat'));rf=load(fullfile(E,'analysis','references.mat'));
tw=rf.tw;old=readtable(fullfile(E,'analysis','scores.csv'),'TextType','string');
runs={'V0','Vm50','Vm30','Vp30','Vp50','R0'};rows={};maxdiff=0;
for i=1:numel(runs)
 L=load(fullfile(E,'results',[runs{i} '.mat']));
 for m={'Input','Standard','QScan_Standard','SC','QSC'}
  if any(strcmp(m{1},{'SC','QSC'}))
   S=load(fullfile(root,'results_sc',sprintf('F15A_%s_%s.mat',runs{i},m{1})));A=S.X;
  else
   A=L.G.(m{1});
  end
  x=path_trace(A,in.idx);t1=tie(zpfilter(x,.004),rf.ref1,tw,rf.w1);t2=tie(x,rf.ref2,tw,rf.w2);
  if ~any(strcmp(m{1},{'SC','QSC'}))
   o=old(old.run==runs{i}&old.method==m{1},:);maxdiff=max([maxdiff abs(o.T1-t1) abs(o.T2-t2)]);
  end
  rows(end+1,:)={runs{i},L.r.e,m{1},t1,t2,L.info.qscan_factor}; %#ok<AGROW>
 end
end
assert(maxdiff<1e-9,'rescoring does not reproduce experiment 5 (max diff %g)',maxdiff);
T=cell2table(rows,'VariableNames',{'run','q_error_percent','method','T1','T2','qscan_factor'});
an=fullfile(root,'analysis');if ~isfolder(an),mkdir(an);end;writetable(T,fullfile(an,'scores_f15a.csv'));
fprintf('experiment-5 scores reproduced (max diff %.2g)\n',maxdiff);disp(T);
end

function x=path_trace(A,idx)
x=zeros(size(A,1),1);for n=1:size(A,1),x(n)=mean(A(n,idx(n)-2:idx(n)+2));end
end
function [best,lag]=tie(x,ref,tw,w)
best=-Inf;lag=NaN;tf=tw(w);rr=ref(w);
for s=-.020:.001:.020
 y=interp1(tw,x,tf+s,'linear',NaN);ok=isfinite(y)&isfinite(rr);a=y(ok)-mean(y(ok));b=rr(ok)-mean(rr(ok));
 c=(a'*b)/sqrt((a'*a)*(b'*b));if c>best,best=c;lag=s*1000;end
end
end
function y=zpfilter(x,dt)
% Experiment 5: |H| = 3rd-order Butterworth high-pass 6 Hz x 12th-order low-pass 70 Hz.
x=x(:);n=numel(x);N=2^nextpow2(2*n);f=abs([0:N/2 -N/2+1:-1]')/(N*dt);
H=1./sqrt(1+(6./max(f,realmin)).^6)./sqrt(1+(f/70).^24);H(1)=0;y=real(ifft(fft(x,N).*H));y=y(1:n);
end
