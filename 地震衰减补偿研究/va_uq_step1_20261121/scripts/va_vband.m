function [vb,chk]=va_vband(M,P,C)
% In-band migration velocity (PLAN_ZH.md section 2): point by point on the true model,
% v_band = v100 [1 + ln(f_band/f_ref)/(pi Q)] with the true Q, then the same Gaussian smoothing
% (sigma P.smooth_m) that make_model uses for M.vmig. chk = max|gsmooth(M.v) - M.vmig| (must be ~0).
fac=1+log(C.f_band/P.fd.fref)./(pi*M.Q);fac(~isfinite(M.Q))=1;
vb=gsmooth(M.v.*fac,P.smooth_m/P.dx);chk=max(abs(reshape(gsmooth(M.v,P.smooth_m/P.dx)-M.vmig,[],1)));
end

function y=gsmooth(x,s)
% copy of make_model.m/gsmooth (local function there)
h=exp(-.5*((-ceil(3*s):ceil(3*s))/s).^2);h=h/sum(h);n=numel(h);m=(n-1)/2;
xp=[repmat(x(1,:),m,1);x;repmat(x(end,:),m,1)];xp=[repmat(xp(:,1),1,m) xp repmat(xp(:,end),1,m)];
y=conv2(h,h,xp,'valid');
end
