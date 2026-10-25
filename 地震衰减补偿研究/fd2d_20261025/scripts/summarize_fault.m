function summarize_fault()
% Near-fault error (PLAN_ZH.md, descriptive item "断层附近的误差"), which the frozen summarize_fd.m
% does not compute. The plan names the item but gives no window, so the windows below were fixed
% on 2026-10-25 before any near-fault number was looked at; nothing in the sweep is rerun.
% Fault plane between traces f and f+1 (fd_geology: x = (f-0.5)*trace_dx, trace j at (j-1)*trace_dx).
%   near(h)  traces f-h+1 .. f+h (h traces on each side), h = 5 primary (±60 m), 3 and 10 as checks
%   far      traces more than 10 from the plane on both sides (j <= f-10 or j >= f+11)
%   all      all 99 traces, a check against analysis/fd_summary.csv
%   cNN      control: every 10-trace block NN..NN+9 lying wholly inside far, to show how much any
%            10-trace window differs from far without a fault (lateral geology, few traces)
% Error in window W: sum_W (X-T)^2 / sum_W T^2, all time samples, averaged over r = 1..6 as in
% summarize_fd.m (normalised squared error). D, B, A as in summarize_fd.m (QScan + Standard,
% nominal Standard, nominal FullCalRefit); noise-free (r = 0) nominal Standard as B0. The
% error-optimal factor is also recomputed inside near(5) and far to see whether it moves.
P=fd_setup();A=fullfile(P.root,'analysis');f=P.factors;nf=numel(f);R=P.R;nx=P.fd.ntraces;
H=[3 5 10];rows={};
for m=1:numel(P.models)
 sp=P.models(m);Tt=load(fullfile(P.root,'truth',[sp.name '_truth.mat']),'target');T=Tt.target;
 j=1:nx;fp=sp.fault_trace;far=j<=fp-10|j>=fp+11;W={true(1,nx),far};wn={'all','far'};  % 'all' must reproduce fd_summary.csv
 for h=H,W{end+1}=j>=fp-h+1&j<=fp+h;wn{end+1}=sprintf('near%d',h);end %#ok<AGROW>
 for j0=1:nx-9,b=j>=j0&j<=j0+9;if all(far(b)),W{end+1}=b;wn{end+1}=sprintf('c%02d',j0);end,end %#ok<AGROW>
 ET=cellfun(@(w) sum(T(:,w).^2,'all'),W);
 for s=1:numel(P.scenes)
  sc=P.scenes(s);dir_=fullfile(P.root,'cache','sweep',[sp.name '_' sc.tag]);
  sq=zeros(R,nf,numel(W));sqA=zeros(R,numel(W));best=zeros(R,1);
  for r=1:R
   L=load(fullfile(dir_,sprintf('r%03d.mat',r)));assert(L.info.optimization_flags==0);best(r)=L.info.best;
   for k=1:nf
    E=(double(L.STD(:,:,k))-T).^2;for w=1:numel(W),sq(r,k,w)=sum(E(:,W{w}),'all')/ET(w);end
   end
   E=(double(L.FCR(:,:,1))-T).^2;for w=1:numel(W),sqA(r,w)=sum(E(:,W{w}),'all')/ET(w);end
  end
  N=load(fullfile(dir_,'r000.mat'));E=(double(N.STD(:,:,1))-T).^2;q=f*sc.q_scale;
  for w=1:numel(W)
   mc=mean(sq(:,:,w),1);[~,io]=min(mc);D=mean(sq(sub2ind([R nf numel(W)],(1:R)',1+best,w*ones(R,1))));
   rows(end+1,:)={sp.name,sc.tag,wn{w},sum(W{w}),D,mc(1),mean(sqA(:,w)),sum(E(:,W{w}),'all')/ET(w),q(io),min(mc)}; %#ok<AGROW>
  end
 end
end
S=cell2table(rows,'VariableNames',{'model','scene','window','n_traces','mse_D','mse_B','mse_A','mse_B_noise_free','q_best','mse_opt'});
writetable(S,fullfile(A,'fault_windows.csv'));disp(S);
end
