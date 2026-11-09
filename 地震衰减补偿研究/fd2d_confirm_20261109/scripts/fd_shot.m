function [att,ela]=fd_shot(Gm,ix0,P,iz_extra)
% Zero-offset shot at grid column ix0 of the global 2-D model Gm: the window of +-half_window
% (plus PML) around it, attenuated and elastic runs, pressure at the source node (and, if
% iz_extra is given, also that many rows straight below it). Ricker source delayed by fd.t0.
fd=P.fd;h=round(fd.half_window/fd.dx)+fd.npml;cols=ix0-h:ix0+h;assert(cols(1)>=1&&cols(end)<=numel(Gm.x));
M=struct('v',Gm.v(:,cols),'rho',Gm.rho(:,cols),'Q',Gm.Q(:,cols));
iz0=fd.npml+1;assert(abs(Gm.z(iz0))<1e-9);src=[iz0 h+1];rec=src;
if nargin>3,rec=[rec;iz0+iz_extra(:) repmat(h+1,numel(iz_extra),1)];end
opt=struct('dx',fd.dx,'dt',fd.dt,'nt',fd.nt,'fref',P.cm.fref,'npml',fd.npml,'src',src,'rec',rec);
att=fd2d_visco(M,ricker_q(P),opt);M.Q(:)=Inf;ela=fd2d_visco(M,ricker_q(P),opt);
end

function q=ricker_q(P)
t=(1:P.fd.nt)'*P.fd.dt-P.fd.t0;a=(pi*P.cm.fpeak*t).^2;q=(1-2*a).*exp(-a);
end
