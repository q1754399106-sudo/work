function [a0,ax,ax_small,p] = ray_absorption(t0,x,c)
% Integral dt/Q along exact two-way ray through horizontal layers.
% The deepest layer extends beneath the last specified interface.
t0=t0(:); n=numel(t0); a0=zeros(n,1); ax=a0; ax_small=a0; p=a0;
bound=[0 cumsum(2*c.h./c.vp(1:3)) Inf];
for k=1:n
    tau=max(0,min(t0(k),bound(2:end))-bound(1:end-1));
    tau(~isfinite(tau))=0;
    active=tau>0;
    if ~any(active)
        % t0=0 is an artificial boundary sample, no reflector is placed here.
        ax(k)=x/(c.vp(1)*c.Q(1)); continue
    end
    v=c.vp(active); q=c.Q(active); tt=tau(active);
    a0(k)=sum(tt./q);
    if x==0, continue; end
    lo=0; hi=(1-1e-12)/max(v);
    for it=1:70
        pm=(lo+hi)/2;
        xx=sum(tt.*v.^2*pm./sqrt(1-(pm*v).^2));
        if xx<x, lo=pm; else, hi=pm; end
    end
    p(k)=(lo+hi)/2;
    ax(k)=sum(tt./q./sqrt(1-(p(k)*v).^2))-a0(k);
    ps=x/sum(tt.*v.^2);
    ax_small(k)=0.5*ps^2*sum(tt.*v.^2./q);
end
end
