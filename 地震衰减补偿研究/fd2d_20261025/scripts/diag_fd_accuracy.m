function diag_fd_accuracy()
% Development diagnostics behind the forward settings (test geology only, none of H1-H6):
%  (1) relaxation frequencies of the three-mechanism GSLS that flatten Q(f) on 5-100 Hz;
%  (2) grid convergence of the elastic zero-offset trace of the test column against kennett_trace,
%      by band, at dx = 5, 4 and 2.5 m.
P=fd_setup();maxNumCompThreads(P.threads);
f=exp(linspace(log(5),log(100),200));w=2*pi*f(:);Qs=[18 40 150];
dev=@(fr) maxdev(fr,w,Qs);
fr0=[1.6 16 160];best=fminsearch(@(u) dev(exp(u)),log([3 20 110]),optimset('Display','off','MaxFunEvals',2000));
fprintf('GSLS f_relax %s: max |Q(f)/Q-1| on 5-100 Hz = %.4f (current %s: %.4f)\n',mat2str(exp(best),4),dev(exp(best)),mat2str(fr0),dev(fr0));
best4=fminsearch(@(u) dev(exp(u)),log([2 9 40 160]),optimset('Display','off','MaxFunEvals',4000));
fprintf('GSLS 4 mechanisms f_relax %s: max dev %.4f\n',mat2str(exp(best4),4),dev(exp(best4)));
for Q=[18 25 40 80 150]
 fprintf('  Q=%d: dev L3 %.4f, L4 %.4f\n',Q,maxdev(exp(best),w,Q),maxdev(exp(best4),w,Q));
end
sp=struct('name','T','seed',2026102799,'macro_vp',[1800 2300 2800 3300],'macro_h',[230 300 350],'macro_Q',[40 60 85 110], ...
 'carbonate',true,'fault_trace',1e6,'throw',0);
c5=[];
for dx=[5 4 2.5]
 Pd=P;Pd.fd.dx=dx;Pd.fd.npml=round(100/dx);Pd.fd.dt=dx*1e-4;Pd.fd.nt=round(1.88/Pd.fd.dt);Pd.fd.half_window=300;
 Gm=fd_geology(sp,Pd);[~,ix50]=min(abs(Gm.x-Gm.trace_x(50)));
 for f_={'v','rho','Q'},Gm.(f_{1})=repmat(Gm.(f_{1})(:,ix50),1,numel(Gm.x));end
 Gm.Q(:)=Inf;H=Gm;H.v(:)=Gm.cap_v;H.rho(:)=Gm.cap_rho;
 [~,e]=fd_shot(Gm,ix50,Pd);[~,h]=fd_shot(H,ix50,Pd,round(Pd.fd.r0/dx));
 fe=process_zo(e,h(:,1),h(:,2),Gm.cap_v,Gm.config,Pd);
 if isempty(c5),c5=Gm.col50;[~,ke]=kennett_trace(c5.v,c5.rho,c5.h,c5.Q,P.cm,P.cm.fpeak);end
 for wdw={[.25 .7],[.7 1.2],[1.2 1.8]}
  k=P.cm.t>=wdw{1}(1)&P.cm.t<wdw{1}(2);x=fe(k);y=ke(k);X=fft(x,512);Y=fft(y,512);fq=(0:256)/(512*.004);X=X(1:257);Y=Y(1:257);
  lo=fq>=10&fq<=40;hi=fq>=50&fq<=80;
  fprintf('dx %.1f m, %.2f-%.2f s: rel %.3f | 10-40 Hz %.3f | 50-80 Hz %.3f\n',dx,wdw{1}(1),wdw{1}(2),norm(x-y)/norm(y), ...
   norm(X(lo)-Y(lo))/norm(Y(lo)),norm(X(hi)-Y(hi))/norm(Y(hi)));
 end
end
end

function d=maxdev(fr,w,Qs)
ts=1./(2*pi*fr(:)');F=sum(1i*w.*ts./(1+1i*w.*ts),2);d=0;
for Q=Qs
 t=fminbnd(@(t) max(abs(qof(t,F)/Q-1)),0,1);d=max(d,max(abs(qof(t,F)/Q-1)));
end
end

function q=qof(t,F)
m=1+t*F;q=real(m)./imag(m);
end
