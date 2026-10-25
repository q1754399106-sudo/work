function [tr,elastic]=kennett_variant(v,rho,h,Q,cm,fpeak,mode)
% kennett_trace with the interaction terms switchable (diagnosis of forward/dictionary mismatch):
%  'full'      exact recursion, all internal multiples and transmission losses (= kennett_trace)
%  'primaries' primaries only, two-way transmission loss (1-r^2) kept
%  'noloss'    primaries only, no transmission loss (convolutional reflectivity)
% Attenuation and dispersion (Kjartansson constant Q) are identical in every mode.
N=cm.nfft;f=(0:N/2)'/(N*cm.dt);w=2*pi*f;w0=2*pi*cm.fref;
src=(2/sqrt(pi))*(f.^2/fpeak^3).*exp(-(f/fpeak).^2);
tr=response(Q);elastic=response(Inf(size(Q)));
 function x=response(q)
  L=numel(v);Rd=zeros(numel(f),1);
  for k=L:-1:2
   zk=impedance(k,q);zu=impedance(k-1,q);r=(zk-zu)./(zk+zu);
   if k<L
    ph=exp(-2i*w.*h(k)./velocity(k,q));
    switch mode
     case 'full',Rd=(r+Rd.*ph)./(1+r.*Rd.*ph);
     case 'primaries',Rd=r+(1-r.^2).*Rd.*ph;
     case 'noloss',Rd=r+Rd.*ph;
     otherwise,error('unknown mode %s',mode);
    end
   else
    Rd=r;
   end
  end
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
