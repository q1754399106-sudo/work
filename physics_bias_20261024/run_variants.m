function run_variants(sets,which)
% which: optional subset of variant rows (post-hoc M5 is run on its own).
% Noise-free QScan selection vs error-optimal factor under forward variants (inputs/PLAN_ZH.md).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(2);
grid=exp(linspace(log(.6),log(2.2),25));traces=1:3:99;rows={};started=tic;
variants={'M1_full','full',false;'M2_primaries','primaries',false;'M3_noloss','noloss',false;'M4_noloss_macroQ','noloss',true;'M5_blocky_sparse','noloss','blocky'};
if nargin<2,which=1:4;end;variants=variants(which,:);tagname=strjoin(variants(:,1)','_');
for s=sets
 in=load(fullfile(root,'dev',[s{1} '.mat']));T=load(fullfile(root,'dev',[s{1} '_truth.mat']),'base','clean');
 cm=in.config;cm.t=cm.t(:);offset=cm.solver_offset;sig=reshape(in.obs(1).sigma(traces),1,[]);base=T.base;nx=99;
 for vv=1:size(variants,1)
  D0=zeros(cm.nt,numel(traces));Tg=D0;
  for q=1:numel(traces)
   [v,rho,h,Qf]=lateral(base,traces(q),nx);
   if isequal(variants{vv,3},'blocky'),[v,rho,h,Qf]=blocky(base,traces(q),nx);
   elseif variants{vv,3},Qf=macro_q(base,traces(q),nx);end
   [D0(:,q),Tg(:,q)]=kennett_variant(v,rho,h,Qf,cm,cm.fpeak,variants{vv,2});
  end
  if strcmp(variants{vv,1},'M1_full')
   rel=norm(D0-T.clean(:,traces),'fro')/norm(T.clean(:,traces),'fro');fprintf('%s M1 reproduction of archived clean: %.3g\n',s{1},rel);
   assert(rel<1e-10,'M1 does not reproduce the archived forward.');
  end
  J=zeros(size(grid));E=J;
  for g=1:numel(grid)
   c=cm;c.Q=cm.Q*grid(g);d=prune_dictionary(c,offset);
   X=refit_lasso(d.Gram,d.B'*D0,.5*sqrt(2*log(cm.nt))*sig,d.L);res=D0-d.B*X;lam=.5*sqrt(2*log(cm.nt))*sig;
   J(g)=sum((.5*sum(res.^2,1)+lam.*sum(abs(X),1))./sig.^2);
   S=zeros(size(D0));for b=1:11:numel(traces),ids=b:min(b+10,numel(traces));o=reflectivity_refit_prepared(D0(:,ids),repmat({d},1,numel(ids)),sig(ids));S(:,ids)=o.gather;end
   E(g)=norm(S-Tg,'fro')/norm(Tg,'fro');
  end
  [~,is]=min(J);[~,io]=min(E);
  rows{end+1}=struct('set',s{1},'variant',variants{vv,1},'c_sel',grid(is),'c_opt',grid(io),'bias_log',log(grid(is)/grid(io)), ...
   'error_at_1',interp1(log(grid),E,0),'error_at_sel',E(is),'error_at_opt',E(io),'grid',grid,'J',J,'E',E); %#ok<AGROW>
  fprintf('%s %-18s c_sel %.3f c_opt %.3f bias %.3f | err(1) %.4f err(sel) %.4f err(opt) %.4f | %.0fs\n',s{1},variants{vv,1},grid(is),grid(io), ...
   log(grid(is)/grid(io)),interp1(log(grid),E,0),E(is),E(io),toc(started));
 end
end
name=sprintf('variants_%s.json',strjoin(sets,'_'));if nargin>=2,name=sprintf('variants_%s_%s.json',strjoin(sets,'_'),tagname);end
uq_writejson(fullfile(root,'analysis',name),[rows{:}]);
end

function [v,rho,h,Q]=blocky(base,j,nx)
% Post-hoc M5: same macro units, thicknesses and travel time, but each unit is
% one homogeneous layer (travel-time average velocity, harmonic Q): sparse reflectivity.
[vf,~,hf,Qf]=lateral(base,j,nx);x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];
hs=base.h.*stretch(base.macro);hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);m=base.macro(hs>0);
hh=hf(2:end);vv=vf(2:numel(hf));qq=Qf(2:numel(hf));v=vf(1);h=hf(1);Q=Qf(1);
for u=unique(m)
 k=m==u;tt=sum(hh(k)./vv(k));v(end+1)=sum(hh(k))/tt;h(end+1)=sum(hh(k));Q(end+1)=tt/sum(hh(k)./vv(k)./qq(k)); %#ok<AGROW>
end
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;
end

function [v,rho,h,Q]=lateral(base,j,nx)
% Identical to make_independent / make_F34 lateral().
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);
end

function Q=macro_q(base,j,nx)
% Each fine layer takes the travel-time harmonic mean Q of its macro unit (smooth Q, same t* per unit).
[v,~,h,Qf]=lateral(base,j,nx);x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];
hs=base.h.*stretch(base.macro);hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);m=base.macro(hs>0);
Q=Qf;tt=h(2:end)./v(2:numel(h));q=Qf(2:numel(h));
for u=unique(m),k=m==u;Q(1+find(k))=sum(tt(k))/sum(tt(k)./q(k));end
Q(end)=Q(end-1);
end
