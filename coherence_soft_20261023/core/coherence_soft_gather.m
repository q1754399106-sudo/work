function out=coherence_soft_gather(R0,supports,dicts,E,V,live,opts)
% Neighbour-coherence-guided soft group shrinkage on fixed FullCalRefit groups.
% gain=E/(E+lambda(c)*V), lambda(c)=lambda_min+(lambda_max-lambda_min)*(1-c)^2.
% R0/supports: unshrunk refit reflectivity and supports (whole gather, nt-by-nx).
% E/V: per-trace group source energy and propagated white-noise energy from
% reflectivity_soft_prepared (unchanged). c is computed from unshrunk group
% contributions of neighbouring traces, never from shrunk output, truth or labels.
% Neighbours only decide the gain; amplitudes stay those of the trace's own refit.
[n,nx]=size(R0);live=logical(reshape(live,1,[]));
gains=cell(1,nx);cs=gains;lambdas=gains;fallback=gains;
contrib=cell(1,nx);centers=cell(1,nx);groups=cell(1,nx);
for j=1:nx
 [groups{j},contrib{j},centers{j}]=group_contributions(R0(:,j),supports{j},dicts{j});
 assert(numel(E{j})==numel(centers{j})&&numel(V{j})==numel(centers{j}),'Group layout differs from soft rule.');
end
R=zeros(n,nx);S=R;
for j=1:nx
 d=dicts{j};ng=numel(centers{j});support=supports{j}(:);
 c=zeros(ng,1);lam=ones(ng,1);fb=true(ng,1);
 near=j+(-opts.half_width:opts.half_width);near(near==j)=[];
 inside=near(near>=1&near<=nx);
 usable=live(j)&&all(live(inside))&&numel(inside)>=2;
 if usable && ng>0
  maxlag=max(0,ceil(opts.dip_periods/(d.frequency*d.dt)));slack=max(1,round(.5/(d.frequency*d.dt)));
  for k=1:ng
   [c(k),~]=group_support(j,k,contrib,centers,inside,maxlag,slack);
   lam(k)=opts.lambda_min+(opts.lambda_max-opts.lambda_min)*(1-c(k))^2;fb(k)=false;
  end
 end
 e=E{j}(:);v=V{j}(:);g=ones(ng,1);
 positive=v>0;g(positive)=e(positive)./(e(positive)+lam(positive).*v(positive));
 assert(all(isfinite(g))&&all(g>=0)&&all(g<=1));
 a=R0(support,j);
 for k=1:ng,a(groups{j}==k)=g(k)*a(groups{j}==k);end
 r=zeros(n,1);r(support)=a;R(:,j)=r;S(:,j)=d.W*r;
 gains{j}=g;cs{j}=c;lambdas{j}=lam;fallback{j}=fb;
end
out=struct('gather',S,'reflectivity',R,'group_gain',{gains},'coherence',{cs}, ...
 'lambda',{lambdas},'fallback',{fallback},'options',opts,'selected_method','CoherenceSoftRefit');
end

function [group,contrib,centers]=group_contributions(r,support,d)
% Same 0.5-period gap rule as reflectivity_soft_prepared.
support=support(:);
if isempty(support),group=[];contrib={};centers=[];return;end
gap=max(1,round(.5/(d.frequency*d.dt)));group=cumsum([1;diff(support)>gap]);ng=max(group);
contrib=cell(ng,1);centers=zeros(ng,1);
for k=1:ng
 ii=support(group==k);w=abs(r(ii));
 contrib{k}=d.W(:,ii)*r(ii);
 centers(k)=sum(ii.*w)/max(sum(w),realmin);
 if sum(w)==0,centers(k)=mean(ii);end
end
end

function [c,detail]=group_support(j,k,contrib,centers,inside,maxlag,slack)
% Per neighbour: best lag-aligned match among the neighbour's groups whose
% centre is inside the dip bound. score=max(0,rho)*min(beta,1/beta), with rho
% the normalised correlation and beta>0 the amplitude projection ratio.
% Side support is the mean over that side; c is the better side, so an event
% that terminates at a fault or pinch-out is still supported from one side.
x=contrib{j}{k};nx=x'*x;score=zeros(1,numel(inside));
if nx<=0,c=0;detail=score;return;end
for q=1:numel(inside)
 m=inside(q);dist=abs(m-j);bound=maxlag*dist;best=0;
 for h=1:numel(centers{m})
  if abs(centers{m}(h)-centers{j}(k))>bound+slack,continue;end
  y=contrib{m}{h};ny=y'*y;if ny<=0,continue;end
  for lag=-bound:bound
   ys=shift(y,lag);xy=x'*ys;if xy<=0,continue;end
   rho=xy/sqrt(nx*ny);beta=xy/nx;
   best=max(best,rho*min(beta,1/beta));
  end
 end
 score(q)=best;
end
sides={inside<j,inside>j};side=zeros(1,2);
for s=1:2
 if any(sides{s}),side(s)=mean(score(sides{s}));end
end
c=min(1,max(side));detail=score;
end

function y=shift(y,lag)
% y(t-lag): positive lag delays the neighbour (zero fill, no wrap).
if lag>0,y=[zeros(lag,1);y(1:end-lag)];elseif lag<0,y=[y(1-lag:end);zeros(-lag,1)];end
end
