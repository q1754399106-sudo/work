function [tr,elastic]=kennett_trace(v,rho,h,Q,cm,fpeak)
% Exact 1D normal-incidence reflection response of a fine-layered viscoacoustic
% stack (no free surface; all internal multiples and transmission losses).
% Kjartansson constant-Q: v(w)=v0*(i w/w0)^gamma, gamma=atan(1/Q)/pi, w0=2*pi*fref.
% Layers 1..L-1 have thickness h; layer L is the lower half-space. The datum
% sits inside layer 1 (no reflection above it). Interfaces are off the time grid.
% Returns attenuated trace and its Q=Inf counterpart, both Ricker-convolved.
N=cm.nfft;f=(0:N/2)'/(N*cm.dt);w=2*pi*f;w0=2*pi*cm.fref;
src=(2/sqrt(pi))*(f.^2/fpeak^3).*exp(-(f/fpeak).^2);   % Ricker amplitude spectrum (zero phase)
tr=response(Q);elastic=response(Inf(size(Q)));
 function x=response(q)
  L=numel(v);Rd=zeros(numel(f),1);
  for k=L:-1:2
   zk=impedance(k,q);zu=impedance(k-1,q);r=(zk-zu)./(zk+zu);
   if k<L
    ph=exp(-2i*w.*h(k)./velocity(k,q));Rd=(r+Rd.*ph)./(1+r.*Rd.*ph);
   else
    Rd=r;
   end
  end
  % Source and receiver at the datum, h(1) above the first interface.
  Rd=Rd.*exp(-2i*w.*h(1)./velocity(1,q));
  X=src.*Rd;X(1)=0;X(end)=real(X(end));X=[X;conj(X(end-1:-1:2))];
  x=real(ifft(X))/cm.dt;x=x(1:cm.nt);
 end
 function c=velocity(k,q)
  if isinf(q(k)),c=v(k)*ones(size(w));return;end
  g=atan(1/q(k))/pi;c=v(k)*(1i*max(w,realmin)/w0).^g;c(1)=v(k);
 end
 function z=impedance(k,q),z=rho(k)*velocity(k,q);end
end
