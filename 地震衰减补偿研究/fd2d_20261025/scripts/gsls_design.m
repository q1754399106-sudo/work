function G=gsls_design(Q,fref)
% Generalised standard linear solid (tau method, Blanch et al. 1995) with four relaxation
% mechanisms, tau chosen per Q to minimise the largest deviation of Q(f) from Q on 5-100 Hz
% (relaxation frequencies from diag_fd_accuracy: max deviation 3.6% at Q = 18, 2.4-3.3% for
% Q = 25-150). Phase velocity is referenced to fref, the reference frequency of the solver
% dictionaries (Kjartansson model, cm.fref = 100 Hz).
%   M(w) = M_R * m(w),  m(w) = 1 + sum_l tau*i*w*ts_l/(1+i*w*ts_l)
% G.ts      relaxation times tau_sigma_l (s), fixed
% G.tau     tau for each entry of Q (0 where Q is Inf)
% G.vfac    Re(1/sqrt(m(w_ref)))^2, so that M_R = rho*v^2*G.vfac gives phase velocity v at fref
% G.q       achieved Q(f) on G.f for each Q (rows), for reporting
G.f_relax=[1.981 8.613 52.86 176.3];G.ts=1./(2*pi*G.f_relax);G.band=[5 100];
G.f=exp(linspace(log(G.band(1)),log(G.band(2)),200));w=2*pi*G.f(:);
F=sum(1i*w.*G.ts./(1+1i*w.*G.ts),2);qof=@(t) real(1+t*F)./imag(1+t*F);
Q=Q(:);tau=zeros(size(Q));G.q=nan(numel(Q),numel(G.f));
for i=find(isfinite(Q))'
 tau(i)=fminbnd(@(t) max(abs(qof(t)/Q(i)-1)),0,2,optimset('TolX',1e-10));G.q(i,:)=qof(tau(i));
end
wr=2*pi*fref;mr=1+tau.*sum(1i*wr*G.ts./(1+1i*wr*G.ts),2);
G.tau=tau;G.vfac=real(1./sqrt(mr)).^2;
end
