function M=make_model(name,P)
% Layered 2-D geology with fine layering and a low-Q lens. Fields (full grid incl. absorbing strips):
%   v     phase velocity at fref (m/s)           Q   true Q (fine layers share their macro layer's Q)
%   vmig  migration velocity (Gaussian-smoothed v) Qm  compensation Q model (macro Q, before any error)
%   rho   constant density
switch name
 case 'P1'   % stage-1 pilot model (development)
  S=struct('seed',2026111701,'z_top',120,'v_top',1850,'Q_top',80, ...
   'iface',{{@(x) 420+.06*(x-1200),@(x) 800-.04*(x-1200),@(x) 1150+.03*(x-1200)}}, ...
   'vp',[2100 2500 2900 3300],'Qp',[55 40 70 100],'lens',[900 600 350 60 18 -.05]);  % xc zc a b Q dv
 case 'T0'   % homogeneous top layer only (direct-wave removal and tests)
  S=struct('seed',1,'z_top',1e9,'v_top',1850,'Q_top',80,'iface',{{}},'vp',[],'Qp',[],'lens',[]);
 otherwise
  error('unknown model %s',name);
end
rng(S.seed,'twister');
[X,Z]=meshgrid(P.x,P.z);Zc=max(Z,0);                         % the top strip copies z = 0
v=S.v_top*ones(size(X));Q=S.Q_top*ones(size(X));
tops=[{@(x) S.z_top+0*x} S.iface];
for k=1:numel(S.vp)
 top=tops{k}(X);
 if k<numel(tops),bot=tops{k+1}(X);else,bot=inf(size(X));end
 in=Zc>=top&Zc<bot;
 % fine layers parallel to the top interface: thickness U(8,40) m, dv N(0,4%) clipped at 8%
 th=8+32*rand(200,1);dv=max(min(.04*randn(200,1),.08),-.08);edges=[0;cumsum(th)];
 d=Zc-top;idx=discretize(d,edges);idx(isnan(idx))=numel(th);
 v(in)=S.vp(k)*(1+dv(idx(in)));Q(in)=S.Qp(k);
end
Qm=Q;
if ~isempty(S.lens)
 L=S.lens;e=((X-L(1))/L(3)).^2+((Zc-L(2))/L(4)).^2<=1;
 v(e)=v(e)*(1+L(6));Q(e)=L(5);Qm(e)=L(5);
end
M=struct('name',name,'v',v,'Q',Q,'Qm',Qm,'rho',2000*ones(size(v)),'vmig',gsmooth(v,P.smooth_m/P.dx),'spec',S);
end

function y=gsmooth(x,s)
h=exp(-.5*((-ceil(3*s):ceil(3*s))/s).^2);h=h/sum(h);n=numel(h);m=(n-1)/2;
xp=[repmat(x(1,:),m,1);x;repmat(x(end,:),m,1)];xp=[repmat(xp(:,1),1,m) xp repmat(xp(:,end),1,m)];
y=conv2(h,h,xp,'valid');
end
