function fe2_smoke(stage)
% Smoke checks (PLAN_ZH.md section 7), no score of any new method.
%  'tab'   : VA operator at Delta t* = 0 is 1; va_tab with known f_v equals uq_tab (no wave solves)
%  'tstar' : anchored Delta t* fields (eikonal only): Delta t* above the anchor, phase of the VA operator
%            just below the anchor, Delta t*/t* in the reservoir window and at F-5
%  'eval'  : fld_eval through fe2_eval on a copy of the stage-4 image file (must reproduce
%            stage4_field.json), fe2_boot on it, fe2_judge on a synthetic json (deleted afterwards)
%  'mig'   : fe2_run('smoke'): one frequency (25 Hz), time and memory
if nargin<1,stage='tab';end
[C,P]=fe2_cfg();
switch stage
 case 'tab'
  xg=single(linspace(0,.05,1201));s=gain_s(70);
  for f=[4 10 25]
   T=va_tab(f,xg,C.sig,C.fv);
   for iv=[0 .1 10]
    H=(T.A(1)+s*T.P(1)*T.C(1)+T.P(1)*iv)/(1+s*T.C(1)+iv);
    fprintf('Delta t* = 0, f %2d Hz, 1/snr %4.1f: H = %.6f%+.6fi\n',f,iv,real(H),imag(H));
   end
   U=uq_tab(f,xg,C.sig,25);V=va_tab(f,xg,C.sig,25);
   fprintf('f %2d: va_tab(f_v = 25 known) vs uq_tab(25): max rel diff A %.2g P %.2g C %.2g\n',f, ...
    max(abs(V.A-U.A))/max(abs(U.A)),max(abs(V.P-U.P))/max(abs(U.P)),max(abs(V.C-U.C))/max(abs(U.C)));
  end
  Lam=log(C.fv(2)/C.fv(1));fprintf('f_v band [%g %g] Hz: geometric centre %.1f Hz, Lambda %.3f; first sinc zero at Delta t* = 0.01 s: %.0f Hz\n',C.fv,sqrt(prod(C.fv)),Lam,pi/(.01*Lam));
 case 'tstar'
  [Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([P.nz P.nx],Zi(:),Xi(:));tic;TS=fe2_tstar(C,P,lin);
  fprintf('t* fields: %.0f s\n',toc);
  z=P.z(Zi(:));z=z(:);x=P.x(Xi(:));x=x(:);za=TS.za_full(Xi(:));za=za(:);ms=1e3;
  ab=z<za;d=double(TS.dlt(ab,:))*ms;
  fprintf('above anchor (%d nodes x %d shots): Delta t* median %.3g ms, 95%% %.3g ms, max %.3g ms, share > 0.1 ms %.3f\n',nnz(ab),size(d,2),median(d(:)),prctile(d(:),95),max(d(:)),mean(d(:)>.1));
  wx=x>=-2500&x<=1500;ab2=ab&wx&z>=za-200;d=double(TS.dlt(ab2,:))*ms;
  fprintf('  within 200 m above the anchor, x -2.5..1.5 km: median %.3g ms, 95%% %.3g ms, max %.3g ms\n',median(d(:)),prctile(d(:),95),max(d(:)));
  r1=z>=za&z<za+10;d=double(TS.dlt(r1,:));
  fprintf('first node row below the anchor: Delta t* median %.3g ms, 95%% %.3g ms, max %.3g ms\n',median(d(:))*ms,prctile(d(:),95)*ms,max(d(:))*ms);
  xg=single(linspace(0,1.05*double(max(TS.abs(:)))/min(C.ms),1201));dxg=double(xg(2));s=gain_s(70);
  for f=[10 25]
   T=va_tab(f,xg,C.sig,C.fv);u=d/min(C.ms)/dxg;i0=min(floor(u),numel(xg)-2);ww=u-i0;i0=i0+1;it=@(t) double(t(i0)).*(1-ww)+double(t(i0+1)).*ww;
   H=(it(T.A)+s*it(T.P).*it(T.C))./(1+s*it(T.C));
   fprintf('  VA operator there (m = 0.5, 1/snr = 0), f %d Hz: max |phase| %.4f rad, |H| %.4f-%.4f\n',f,max(abs(angle(H(:)))),min(abs(H(:))),max(abs(H(:))));
  end
  wr=z>=2600&z<=3200&wx;a=double(TS.abs(wr,:));dd=double(TS.dlt(wr,:));q=dd(:)./a(:);
  fprintf('reservoir window: t* median %.2f ms, Delta t* median %.2f ms, Delta t*/t* median %.3f (5-95%%: %.3f-%.3f)\n',median(a(:))*ms,median(dd(:))*ms,median(q),prctile(q,5),prctile(q,95));
  wf=x>=-1280&x<=-1110&z>=2570&z<=3170;a=double(TS.abs(wf,:));dd=double(TS.dlt(wf,:));
  fprintf('F-5 tie columns 2570-3170 m: t* median %.2f ms, Delta t* median %.2f ms, ratio %.3f\n',median(a(:))*ms,median(dd(:))*ms,median(dd(:)./a(:)));
  fprintf('max Delta t* %.2f ms, max t* %.2f ms\n',max(TS.dlt(:))*ms,max(TS.abs(:))*ms);
 case 'eval'
  copyfile(fullfile(C.s4cache,'images.mat'),fullfile(C.cache,'images_fe2_evaltest.mat'));
  fe2_eval('_evaltest');
  A=jsondecode(fileread(fullfile(C.root,'analysis','fe2_eval_evaltest.json')));
  B=jsondecode(fileread(fullfile(C.s4,'analysis','stage4_field.json')));md=abs(A.R0.tie_F5-B.R0.tie_F5);
  for q={'PH','IQ30','IQ70','LP','W0','UQW0','W1','UQW1'},md=max([md;abs(A.(q{1}).tie_F5-B.(q{1}).tie_F5);abs(A.(q{1}).cen(:)-B.(q{1}).cen(:))]);end
  fprintf('fe2_eval on the stage-4 image file vs stage4_field.json: max |diff| %.3g (R0 tie %.4f)\n',md,A.R0.tie_F5);
  fe2_boot('_evaltest');
  % judge on a synthetic json (VA set to 0.58..0.64, R0 from stage 4), files removed afterwards
  A.VA=struct('tie_F5',linspace(.58,.64,9)');fj=fullfile(C.root,'analysis','fe2_eval_judgetest.json');
  fid=fopen(fj,'w');fprintf(fid,'%s',jsonencode(A));fclose(fid);J=fe2_judge('_judgetest');
  fprintf('judge test (expected P4 0, R0+-0.03 0, follow-up 0): %d %d %d\n',J.P4,J.P4_R0_pm_003,J.followup_2D_prestack_diagnostics);
  delete(fj);delete(fullfile(C.root,'analysis','fe2_judge_judgetest.json'));
 case 'mig'
  fe2_run('smoke');
end
end
