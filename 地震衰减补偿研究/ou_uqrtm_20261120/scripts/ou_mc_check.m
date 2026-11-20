function T=ou_mc_check(name,nmc)
% Step 0: Monte-Carlo check of the path-integral variance (ou_sigma) against simulated fields
% du with separable exponential covariance (AR(1) filters along z and x), and of the simple
% zero-offset formula sigma_eff^2 = (2 Lz / z) s^2 (uniform a, uniform s). Writes analysis/mc_check_<name>.json.
if nargin<1,name='D2';end
if nargin<2,nmc=400;end
P=ou_setup();maxNumCompThreads(P.threads);M=make_model(name,P);ou=P.ou;
Fl=ou_flags(M,P);sloc=ou.s_bg+(ou.s_flag-ou.s_bg)*Fl;
pts=[1000 1200;1000 600;500 1200;1400 1800;800 1600];               % z, x (m)
offs=[0 800];                                                       % |x - x_s| (m)
xs=P.x(P.ixs);zs=P.z(P.izs);
cases=[];for i=1:size(pts,1),for o=offs,[~,k]=min(abs(xs-(pts(i,2)-o)));cases(end+1,:)=[pts(i,:) k]; end,end %#ok<AGROW>
% formula on the same points: ou_sigma on a tiny grid is not possible (it maps the interior), so
% evaluate the leg sums directly with the same code path through a 1-shot wrapper
Rf=ou_sigma(P,M.Qm,M.vmig,sloc,ou);zi=P.z(P.in{1});xi=P.x(P.in{2});
rf=zeros(size(cases,1),1);
for c=1:size(cases,1)
 [~,iz]=min(abs(zi-cases(c,1)));[~,ix]=min(abs(xi-cases(c,2)));rf(c)=Rf.r(sub2ind([numel(zi) numel(xi)],iz,ix),cases(c,3));
end
% Monte Carlo: fine sampling of the legs
a=1./(M.Qm.*M.vmig);Ga=griddedInterpolant({P.z,P.x},a,'linear','nearest');
nf=600;tau=((1:nf)-.5)/nf;rng(20261120,'twister');
pz=exp(-P.dx/ou.Lz);px=exp(-P.dx/ou.Lx);tt=zeros(size(cases,1),nmc);
legs=cell(size(cases,1),1);
for c=1:size(cases,1)
 z=cases(c,1);x=cases(c,2);k=cases(c,3);xr=min(max(2*x-xs(k),xs(1)),xs(end));
 q=[z+(zs-z)*tau z+(zs-z)*tau;x+(xs(k)-x)*tau x+(xr-x)*tau];
 ds=[hypot(z-zs,x-xs(k))/nf*ones(1,nf) hypot(z-zs,x-xr)/nf*ones(1,nf)];
 legs{c}=struct('q',q,'w',Ga(q(1,:),q(2,:)).*ds);
end
for r=1:nmc
 e=randn(P.nz,P.nx);
 for i=2:P.nz,e(i,:)=pz*e(i-1,:)+sqrt(1-pz^2)*e(i,:);end
 for j=2:P.nx,e(:,j)=px*e(:,j-1)+sqrt(1-px^2)*e(:,j);end
 du=sloc.*e;Gd=griddedInterpolant({P.z,P.x},du,'linear','nearest');
 for c=1:size(cases,1),L=legs{c};tt(c,r)=sum(L.w.*(1+Gd(L.q(1,:),L.q(2,:))));end
end
t0=zeros(size(cases,1),1);for c=1:size(cases,1),t0(c)=sum(legs{c}.w);end
rmc=var(tt,0,2)./t0.^2;
% simple formula: uniform s = s_bg, zero offset, 2 Lz / z
rs=2*ou.Lz./cases(:,1)*ou.s_bg^2;
T=struct('model',name,'nmc',nmc,'z',cases(:,1),'x',cases(:,2),'xs',xs(cases(:,3))','r_formula',rf,'r_mc',rmc, ...
 'r_mc_se',rmc*sqrt(2/(nmc-1)),'sig_formula',sqrt(log1p(rf)),'sig_mc',sqrt(log1p(rmc)),'sig_simple_bg',sqrt(rs),'ou',ou);
disp(table(T.z,T.x,T.xs,T.r_formula,T.r_mc,T.r_mc_se,T.sig_formula,T.sig_mc,T.sig_simple_bg));
fid=fopen(fullfile(P.root,'analysis',['mc_check_' name '.json']),'w');fprintf(fid,'%s',jsonencode(T,'PrettyPrint',true));fclose(fid);
end
