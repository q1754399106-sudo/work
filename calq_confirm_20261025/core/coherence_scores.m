function sc=coherence_scores(R0,supports,d,live,half_width,dip_periods)
% Neighbour support c per group, identical to coherence_soft_gather (v1), plus
% the best-match lag to each neighbour (for later stages). One dictionary d
% (single section offset). Unshrunk contributions only; no truth.
[~,nx]=size(R0);live=logical(reshape(live,1,[]));
gap=max(1,round(.5/(d.frequency*d.dt)));maxlag=max(0,ceil(dip_periods/(d.frequency*d.dt)));
groups=cell(1,nx);contrib=groups;centers=groups;
for j=1:nx
 s=supports{j}(:);if isempty(s),groups{j}=[];contrib{j}={};centers{j}=zeros(0,1);continue;end
 g=cumsum([1;diff(s)>gap]);ng=max(g);contrib{j}=cell(ng,1);centers{j}=zeros(ng,1);groups{j}=g;
 for k=1:ng
  ii=s(g==k);w=abs(R0(ii,j));contrib{j}{k}=d.W(:,ii)*R0(ii,j);
  centers{j}(k)=sum(ii.*w)/max(sum(w),realmin);if sum(w)==0,centers{j}(k)=mean(ii);end
 end
end
c=cell(1,nx);fallback=c;lags=c;scores=c;
for j=1:nx
 ng=numel(centers{j});c{j}=zeros(ng,1);fallback{j}=true(ng,1);
 near=j+(-half_width:half_width);near(near==j)=[];inside=near(near>=1&near<=nx);
 lags{j}=nan(ng,numel(near));scores{j}=zeros(ng,numel(near));
 if ~(live(j)&&all(live(inside))&&numel(inside)>=2)||ng==0,continue;end
 for k=1:ng
  x=contrib{j}{k};nxx=x'*x;fallback{j}(k)=false;if nxx<=0,continue;end
  sc_=zeros(1,numel(inside));lg=nan(1,numel(inside));
  for q=1:numel(inside)
   m=inside(q);bound=maxlag*abs(m-j);best=0;bl=NaN;
   for h=1:numel(centers{m})
    if abs(centers{m}(h)-centers{j}(k))>bound+gap,continue;end
    y=contrib{m}{h};ny=y'*y;if ny<=0,continue;end
    for lag=-bound:bound
     ys=shiftv(y,lag);xy=x'*ys;if xy<=0,continue;end
     v=(xy/sqrt(nxx*ny))*min(xy/nxx,nxx/xy);
     if v>best,best=v;bl=lag+(centers{m}(h)-centers{j}(k));end
    end
   end
   sc_(q)=best;lg(q)=bl;
  end
  side=[0 0];L=inside<j;Rr=inside>j;
  if any(L),side(1)=mean(sc_(L));end
  if any(Rr),side(2)=mean(sc_(Rr));end
  c{j}(k)=min(1,max(side));
  pos=ismember(near,inside);scores{j}(k,pos)=sc_;lags{j}(k,pos)=lg;
 end
end
sc=struct('c',{c},'fallback',{fallback},'groups',{groups},'centers',{centers},'lags',{lags},'scores',{scores},'gap',gap);
end

function y=shiftv(y,lag)
if lag>0,y=[zeros(lag,1);y(1:end-lag)];elseif lag<0,y=[y(1-lag:end);zeros(-lag,1)];end
end
