function [fault,pinch,lowE,deep]=complex_regions(model,g)
% Truth-derived scoring masks (evaluation only); same definitions as v1 V3.
t=model.t;x=model.x;[nt,nx]=size(model.datasets(1).target);mt=model.truth;
fault=false(nt,nx);pinch=false(nt,nx);deep=false(nt,nx);deep(t>=1.05,model.live)=true;
for j=1:nx
 if any(abs(x(j)-g.faults)<=90),fault(t>=mt.cumtime(61,j),j)=true;end
 if x(j)>=g.wedge_tip-450&&x(j)<=g.wedge_tip+450
  iz=find(mt.strata(:,j)>=g.wedge_depth,1);if ~isempty(iz),pinch(abs(t-mt.cumtime(iz,j))<=.04,j)=true;end
 end
end
fault(:,~model.live)=false;pinch(:,~model.live)=false;lowE=false(nt,nx);
for k=1:numel(model.datasets)
 T=model.datasets(k).target;L=sqrt(movmean(T.^2,21,1));rowrms=sqrt(mean(T.^2,2));sig=rowrms>.01*max(rowrms);
 vals=L(sig,model.live);lo=prctile(vals(:),20);hi=prctile(vals(:),45);
 m=false(nt,nx);m(sig,model.live)=L(sig,model.live)>=lo&L(sig,model.live)<=hi;lowE(:,:,k)=m;
end
end
