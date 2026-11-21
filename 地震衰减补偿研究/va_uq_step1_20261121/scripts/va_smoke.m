function va_smoke(stage)
% Smoke checks (PLAN_ZH.md section 7). stage 'tab': tables, sinc factor, GSLS dispersion, v_band (no
% wave solves); 'tstar': t* fields and anchor (eikonal only); 'mig': one frequency of va_mig on D2.
if nargin<1,stage='tab';end
[C,P]=va_cfg('D2');
switch stage
 case 'tab'
  xg=single(linspace(0,.05,1201));
  for f=[10 25 40 55]
   T=va_tab(f,xg,C.sig,C.fv);s=gain_s(70);
   for iv=[0 .1 10]
    H=(T.A(1)+s*T.P(1)*T.C(1)+T.P(1)*iv)/(1+s*T.C(1)+iv);
    fprintf('anchor (Delta t* = 0) f %2d Hz 1/snr %4.1f: H = %.6f%+.6fi\n',f,iv,real(H),imag(H));
   end
   U=uq_tab(f,xg,C.sig,C.f_band);V=va_tab(f,xg,C.sig,C.f_band);
   fprintf('f %2d: va_tab(fixed 25 Hz) vs uq_tab(25 Hz): max rel diff A %.2g P %.2g C %.2g\n',f, ...
    max(abs(V.A-U.A))/max(abs(U.A)),max(abs(V.P-U.P))/max(abs(U.P)),max(abs(V.C-U.C))/max(abs(U.C)));
  end
  % closed-form f_v expectation against brute force (sig = 0, so x is exact)
  lf=linspace(log(C.fv(1)),log(C.fv(2)),20001);
  for f=[10 25 55]
   for x=[.002 .01 .05]
    bf=mean(exp(-1i*2*f*x*(log(f)-lf)));T=va_tab(f,x,0,C.fv);
    fprintf('E_fv phase f %2d x %.3f: closed %.6f%+.6fi brute %.6f%+.6fi sinc %.5f\n',f,x,real(T.P),imag(T.P),real(bf),imag(bf),sin(f*x*C.Lam)/(f*x*C.Lam));
   end
  end
  fprintf('first sinc zero at Delta t* = 0.01 s: f = %.0f Hz\n',pi/(.01*C.Lam));
  % GSLS phase velocity ratio v(25)/v(100) against the constant-Q formula used for v_band
  Ql=[20 22 50 60 70 80 90 110 120];G=gsls_design(Ql,P.fd.fref);
  m=@(f) 1+G.tau.*sum(1i*2*pi*f*G.ts./(1+1i*2*pi*f*G.ts),2);vr=real(1./sqrt(m(P.fd.fref)))./real(1./sqrt(m(C.f_band)));
  for i=1:numel(Ql),fprintf('Q %3d: GSLS v25/v100 = %.5f, formula %.5f\n',Ql(i),vr(i),1+log(C.f_band/P.fd.fref)/(pi*Ql(i)));end
  for n={'D2','D4'}
   M=make_model(n{1},P);[vb,chk]=va_vband(M,P,C);rr=vb./M.vmig;rr=rr(P.in{1},P.in{2});
   fprintf('%s v_band: smoothing check %.3g m/s, vb/vmig interior min %.4f median %.4f max %.4f\n',n{1},chk,min(rr(:)),median(rr(:)),max(rr(:)));
  end
 case 'tstar'
  for n={'D2','D4'}
   [C,P]=va_cfg(n{1});M=make_model(n{1},P);vb=va_vband(M,P,C);
   [Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([P.nz P.nx],Zi(:),Xi(:));tic;TS=va_tstar(M,P,C,vb,lin,Xi(:));
   z=P.z(Zi(:));x=P.x(Xi(:));zaq=TS.za(Xi(:));zaq=zaq(:);
   at=z(:)>=zaq&z(:)<zaq+P.dx&x(:)>=300&x(:)<=2100;dp=z(:)>=820&z(:)<=1550&x(:)>=300&x(:)<=2100;
   ab=z(:)<zaq-P.dx&x(:)>=300&x(:)<=2100;ra=TS.dlt(dp,:)./TS.abs(dp,:);
   da=TS.dlt(at,:);dab=TS.dlt(ab,:);
   fprintf('%s (%.0f s): anchor %.0f-%.0f m; first node row below anchor: Delta t* median %.2g max %.2g s (t* there median %.4f s)\n', ...
    n{1},toc,min(TS.za(P.in{2})),max(TS.za(P.in{2})),median(da(:)),max(da(:)),median(TS.abs(at,:),'all'));
   fprintf('%s above anchor (x 300-2100): Delta t* median %.2g, 99%% %.2g, max %.2g s; fraction > 1e-4 s: %.4f\n',n{1}, ...
    median(dab(:)),prctile(dab(:),99),max(dab(:)),mean(dab(:)>1e-4));
   fprintf('%s W_deep: t* median %.4f s, Delta t* median %.4f s, Delta/t* median %.3f (5-95%%: %.3f-%.3f)\n',n{1}, ...
    median(TS.abs(dp,:),'all'),median(TS.dlt(dp,:),'all'),median(ra,'all'),prctile(ra(:),5),prctile(ra(:),95));
  end
 case 'mig'
  j=find(P.f==25);va_mig('D2',j);
  I=load(fullfile(P.cache,'images_va_D2_smoke.mat'));
  for q=[{'R0'} I.meta.methods],v=I.(q{1});fprintf('%-9s rms %s\n',q{1},mat2str(sqrt(mean(double(v).^2,1)),3));end
end
end
