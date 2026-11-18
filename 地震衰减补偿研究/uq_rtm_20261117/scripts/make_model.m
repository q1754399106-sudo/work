function M=make_model(name,P)
% Layered 2-D geology with fine layering and optional low-Q bodies. Fields (full grid incl. absorbing strips):
%   v     phase velocity at fref (m/s)           Q   true Q (fine layers share their macro layer's Q)
%   vmig  migration velocity (Gaussian-smoothed v) Qm  compensation Q model (macro Q and bodies, before any error)
%   rho   constant density
% Optional features (absent = off; P1 and T0 are unchanged by them, bit for bit):
%   lens    rows [xc zc a b Q dv]: elliptical bodies          fault  [x0 dip_deg throw]: normal fault below z_top,
%   chimney [xc halfwidth ztop zbot Q dv]: vertical column       hanging wall (x > plane) shifted down by throw
%   qgrad   [fl fr]: Q below z_top multiplied by a factor going linearly from fl (x = 0) to fr (x = 2400 m)
switch name
 case 'P1'   % stage-1 pilot model (development)
  S=struct('seed',2026111701,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 420+.06*(x-1200),@(x) 800-.04*(x-1200),@(x) 1150+.03*(x-1200)}}, ...
   'vp',[2100 2500 2900 3300],'Qp',[55 40 70 100],'lens',[900 600 350 60 18 -.05]);  % xc zc a b Q dv
 case 'T0'   % homogeneous top layer only (direct-wave removal and tests)
  S=struct('seed',1,'z_top',1e9,'v_top',1850,'Q_top',80,'iface',{{}},'vp',[],'Qp',[],'lens',[]);
 % ---- stage-2 development models (PLAN_STAGE2_DRAFT_ZH.md) ----
 case 'D1'   % fault and two lenses
  S=struct('seed',2026111711,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 380+.05*(x-1200),@(x) 760+.02*(x-1200),@(x) 1120-.03*(x-1200)}}, ...
   'vp',[2050 2450 2850 3250],'Qp',[60 45 75 110],'lens',[700 520 250 50 20 -.04;1700 900 300 60 25 -.04], ...
   'fault',[1300 70 40]);
 case 'D2'   % anticline and gas chimney
  S=struct('seed',2026111721,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 450-120*exp(-((x-1200)/500).^2),@(x) 820-150*exp(-((x-1200)/550).^2),@(x) 1180-100*exp(-((x-1200)/600).^2)}}, ...
   'vp',[2150 2550 2950 3350],'Qp',[70 50 60 90],'lens',[],'chimney',[1200 120 300 800 22 -.03]);
 case 'D3'   % laterally graded Q
  S=struct('seed',2026111731,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 400+.03*(x-1200),@(x) 780+.03*(x-1200),@(x) 1150+.03*(x-1200)}}, ...
   'vp',[2100 2500 2900 3300],'Qp',[60 45 70 100],'lens',[],'qgrad',[1.6 .6]);
 case 'D4'   % shallow thin low-Q layer and steep dips
  S=struct('seed',2026111741,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 260+0*x,@(x) 330+0*x,@(x) 600+.25*(x-1200),@(x) 1200-.15*(x-1200)}}, ...
   'vp',[2000 2250 2700 3200 3400],'Qp',[20 60 80 110 120],'lens',[]);
 % ---- stage-3 confirmation models (PLAN_STAGE3_ZH.md), never used in development ----
 case 'C1'   % fault and gas chimney
  S=struct('seed',2026111901,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 360+.04*(x-1200),@(x) 740-.03*(x-1200),@(x) 1100+.05*(x-1200)}}, ...
   'vp',[2000 2400 2800 3200],'Qp',[65 50 80 105],'lens',[],'fault',[900 65 35],'chimney',[1600 100 350 750 25 -.03]);
 case 'C2'   % three stacked lenses and a reversed lateral Q gradient
  S=struct('seed',2026111911,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 420+0*x,@(x) 790+.05*(x-1200),@(x) 1160+.02*(x-1200)}}, ...
   'vp',[2150 2500 2950 3350],'Qp',[70 55 65 95], ...
   'lens',[800 550 280 45 22 -.04;1300 900 300 50 20 -.04;1800 650 250 40 28 -.03],'qgrad',[.7 1.4]);
 case 'C3'   % syncline and a small shallow low-Q lens
  S=struct('seed',2026111921,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 400+130*exp(-((x-1200)/550).^2),@(x) 780+160*exp(-((x-1200)/600).^2),@(x) 1150+110*exp(-((x-1200)/650).^2)}}, ...
   'vp',[2050 2450 2850 3250],'Qp',[60 45 70 90],'lens',[1000 300 200 35 15 -.05]);
 case 'C4'   % high-Q overburden over a deep low-Q layer
  S=struct('seed',2026111931,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 380+.04*(x-1200),@(x) 760+.04*(x-1200),@(x) 980+.04*(x-1200),@(x) 1200+.04*(x-1200)}}, ...
   'vp',[2100 2600 2900 3100 3400],'Qp',[90 110 30 85 120],'lens',[]);
 otherwise
  error('unknown model %s',name);
end
rng(S.seed,'twister');
[X,Z]=meshgrid(P.x,P.z);Zc=max(Z,0);                         % the top strip copies z = 0
Zg=Zc;                                                       % geological depth (fault-shifted)
if isfield(S,'fault')&&~isempty(S.fault)
 F=S.fault;hw=X>F(1)+Zc/tan(F(2)*pi/180)&Zc>S.z_top;Zg(hw)=Zc(hw)-F(3);
end
v=S.v_top*ones(size(X));Q=S.Q_top*ones(size(X));
tops=[{@(x) S.z_top+0*x} S.iface];
for k=1:numel(S.vp)
 top=tops{k}(X);
 if k<numel(tops),bot=tops{k+1}(X);else,bot=inf(size(X));end
 in=Zg>=top&Zg<bot;
 % fine layers parallel to the top interface: thickness U(8,40) m, dv N(0,4%) clipped at 8%
 th=8+32*rand(200,1);dv=max(min(.04*randn(200,1),.08),-.08);edges=[0;cumsum(th)];
 d=Zg-top;idx=discretize(d,edges);idx(isnan(idx))=numel(th);
 v(in)=S.vp(k)*(1+dv(idx(in)));Q(in)=S.Qp(k);
end
if isfield(S,'qgrad')&&~isempty(S.qgrad)
 g=S.qgrad(1)+(S.qgrad(2)-S.qgrad(1))*min(max(X/2400,0),1);b=Zc>S.z_top;Q(b)=Q(b).*g(b);
end
Qm=Q;
for i=1:size(S.lens,1)
 L=S.lens(i,:);e=((X-L(1))/L(3)).^2+((Zc-L(2))/L(4)).^2<=1;
 v(e)=v(e)*(1+L(6));Q(e)=L(5);Qm(e)=L(5);
end
if isfield(S,'chimney')&&~isempty(S.chimney)
 C=S.chimney;e=abs(X-C(1))<=C(2)&Zc>=C(3)&Zc<=C(4);v(e)=v(e)*(1+C(6));Q(e)=C(5);Qm(e)=C(5);
end
M=struct('name',name,'v',v,'Q',Q,'Qm',Qm,'rho',2000*ones(size(v)),'vmig',gsmooth(v,P.smooth_m/P.dx),'spec',S);
end

function y=gsmooth(x,s)
h=exp(-.5*((-ceil(3*s):ceil(3*s))/s).^2);h=h/sum(h);n=numel(h);m=(n-1)/2;
xp=[repmat(x(1,:),m,1);x;repmat(x(end,:),m,1)];xp=[repmat(xp(:,1),1,m) xp repmat(xp(:,end),1,m)];
y=conv2(h,h,xp,'valid');
end
