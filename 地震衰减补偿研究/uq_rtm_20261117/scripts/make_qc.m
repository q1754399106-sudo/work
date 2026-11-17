function Qc=make_qc(M,P,form,m)
% Compensation Q model of PLAN_STAGE2_ZH.md section 2 (full grid). For E1-E3 call with m = 1 and
% scale t* by 1/m afterwards; E4 is not a pure scaling and takes its own m.
[X,Z]=meshgrid(P.x,P.z);Zc=max(Z,0);Qm=M.Qm;below=Zc>=M.spec.z_top;
switch form
 case 'E1'
  Qc=Qm*m;
 case 'E2'                                         % mean 1/Q below z_top over the interior
  I=false(size(Qm));I(P.in{1},P.in{2})=true;Qc=Qm;
  Qc(below)=1/mean(1./Qm(below&I));Qc=Qc*m;
 case 'E3'                                         % Gaussian smoothing of 1/Q, sigma 300 m (x), 50 m (z)
  Qc=1./gsm(1./Qm,50/P.dx,300/P.dx)*m;
 case 'E4'                                         % Q m^(2z/1600 - 1)
  Qc=Qm.*m.^(2*min(Zc,1600)/1600-1);
 otherwise
  error('unknown form %s',form);
end
end

function y=gsm(x,sz,sx)
hz=exp(-.5*((-ceil(3*sz):ceil(3*sz))/sz).^2);hz=hz/sum(hz);hx=exp(-.5*((-ceil(3*sx):ceil(3*sx))/sx).^2);hx=hx/sum(hx);
mz=(numel(hz)-1)/2;mx=(numel(hx)-1)/2;
xp=[repmat(x(1,:),mz,1);x;repmat(x(end,:),mz,1)];xp=[repmat(xp(:,1),1,mx) xp repmat(xp(:,end),1,mx)];
y=conv2(hz,hx,xp,'valid');
end
