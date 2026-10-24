function predict()
% Truth-free local output sensitivity diagnostic. No targets, no new solver.
root=fileparts(fileparts(mfilename('fullpath')));repo=fileparts(root);
addpath(fullfile(repo,'bae_q_20261024','core'),fullfile(repo,'bae_q_20261024','source'));
maxNumCompThreads(4);rng(20261024,'twister');started=tic;
files=dir(fullfile(root,'inputs','*.mat'));rows={};checks={};audits={};
for fi=1:numel(files)
 s=load(fullfile(files(fi).folder,files(fi).name));
 assert(~any(isfield(s,{'target','truth','reflectivity','tr'})));
 c=s.cm;off=450;if isfield(c,'solver_offset'),off=c.solver_offset;end
 [K,caldiag]=calibration_information(s.calibration,c);op=operators(c,off,0,true);
 scan=c;scan.Q=s.scaninfo.applied_Q;ops=operators(scan,off,0,true);
 d=prune_dictionary(c,off);ds=prune_dictionary(scan,off);
 parity=norm(op.A-d.B.*d.scale','fro')/norm(op.A,'fro');assert(parity<1e-11);
 cache=fullfile(root,'cache');if ~isfolder(cache),mkdir(cache);end
 cached=fullfile(cache,[s.name '_' s.condition '_supports.mat']);
 if isfile(cached)
  z=load(cached);R=z.R;Rs=z.Rs;support=z.support;supports=z.supports;recomputed=z.recomputed;recomputed_scan=z.recomputed_scan;
 else
  R=zeros(size(s.D));Rs=R;recomputed=R;recomputed_scan=R;support=cell(1,size(R,2));supports=support;
  for first=1:11:size(s.D,2)
   ids=first:min(first+10,size(s.D,2));
   o=reflectivity_refit_prepared(s.D(:,ids),repmat({d},1,numel(ids)),s.sigma(ids));
   os=reflectivity_refit_prepared(s.D(:,ids),repmat({ds},1,numel(ids)),s.sigma(ids));
   R(:,ids)=o.reflectivity;Rs(:,ids)=os.reflectivity;support(ids)=o.support;supports(ids)=os.support;
   recomputed(:,ids)=o.gather;recomputed_scan(:,ids)=os.gather;
  end
  save(cached,'R','Rs','support','supports','recomputed','recomputed_scan','-v7');
 end
 base_diff=norm(recomputed-s.FullCalRefit,'fro')/max(norm(s.FullCalRefit,'fro'),realmin);
 scan_diff=norm(recomputed_scan-s.QScanRefit,'fro')/max(norm(s.QScanRefit,'fro'),realmin);
 assert(max(base_diff,scan_diff)<1e-6,'Frozen algorithm reproduction failed.');
 checks=[checks;derivative_checks(op,c,off,R(:,s.trace_ids(1)),s.name,s.condition)]; %#ok<AGROW>
 [Kscan,~]=calibration_information(s.calibration,scan);
 edges=.30:.05:1.60;
 for j=s.trace_ids
  sig=max(s.sigma(j),realmin);active=any(s.D(:,j)~=0);
  ids=support{j};idscan=supports{j};
  [J,T]=jacobians(op,R(:,j),ids,sig);[Js,Ts]=jacobians(ops,Rs(:,j),idscan,sig);
  k=size(J,2)-6;ks=size(Js,2)-6;
  [ur,ir,~]=risk(J(:,1:k),T(:,1:k),1e-10);
  [uq,iq,~]=risk(J(:,1:k+4),T(:,1:k+4),1e-10);
  [uf,ifull,~]=risk(J,T,1e-10);
  Ja=[J;zeros(size(K,1),k),K];Jsa=[Js;zeros(size(Kscan,1),ks),Kscan];
  [uc,ic,rankc]=risk(Ja,T,1e-10);[us,is,~]=risk(Jsa,Ts,1e-10);
  [ut8,it8,~]=risk(Ja,T,1e-8);[ut12,it12,~]=risk(Ja,T,1e-12);
  dense=union(ids,1:max(1,round(.012/c.dt)):c.nt);
  [Je,Te]=jacobians(op,R(:,j),dense,sig);ke=size(Je,2)-6;
  [ue,ie,ranke]=risk([Je;zeros(size(K,1),ke),K],Te,1e-10);
  audits(end+1,:)={s.name,s.condition,j,active,numel(ids),numel(dense),rankc,ranke,base_diff,scan_diff,parity, ...
   size(K,1),caldiag.reduced_residual,caldiag.rank,mean(ic),mean(ie)}; %#ok<AGROW>
  for w=1:numel(edges)-1
   mask=c.t>=edges(w)-1e-9&c.t<edges(w+1)-1e-9;
   obsrms=sqrt(mean(s.D(mask,j).^2));snrrisk=sig/max(obsrms,sig*1e-6);
   spread=sqrt(mean(s.qspread(mask,j).^2));
   vals={s.name,s.condition,j,w,edges(w),edges(w+1),active,(edges(w)+edges(w+1))/2,snrrisk,spread, ...
    wrms(ur,ir,mask),wrms(uq,iq,mask),wrms(uf,ifull,mask),wrms(uc,ic,mask), ...
    wrms(ue,ie,mask),wrms(us,is,mask),wrms(ut8,it8,mask),wrms(ut12,it12,mask), ...
    sqrt(mean(uc(mask).^2)),sqrt(mean(ue(mask).^2)),any(ic(mask)),any(ie(mask))};
   rows(end+1,:)=vals; %#ok<AGROW>
  end
 end
 fprintf('%s %s: %d traces, baseline diff %.3g, scan diff %.3g, calibration rank %d, elapsed %.1fs\n', ...
  s.name,s.condition,numel(s.trace_ids),base_diff,scan_diff,caldiag.rank,toc(started));
 writetable(cell2table(rows,'VariableNames',{'model','condition','trace','window','start_s','end_s','live','depth_time', ...
  'inverse_observed_snr','q_spread','U_r','U_Q','U_Q_f_phase','U_cal','U_expanded','U_scan','U_cal_tol8','U_cal_tol12', ...
  'finite_part_cal','finite_part_expanded','null_cal','null_expanded'}),fullfile(root,'predictions','risk.csv'));
end
writetable(cell2table(checks,'VariableNames',{'model','condition','parameter','step','observation_relative','target_relative','target_absolute'}), ...
 fullfile(root,'analysis','derivative_checks.csv'));
writetable(cell2table(audits,'VariableNames',{'model','condition','trace','live','support_size','expanded_size','joint_rank','expanded_rank', ...
 'frozen_base_difference','frozen_scan_difference','operator_parity','calibration_rows','calibration_reduced_residual', ...
 'calibration_rank','null_fraction','expanded_null_fraction'}),fullfile(root,'analysis','numerical_audit.csv'));
fprintf('DONE prediction %.1fs\n',toc(started));
end

function op=operators(c,offset,phi,derivatives)
t=c.t(:);z=pi*c.fpeak*(t-t');W=(1-2*z.^2).*exp(-z.^2);
Wp=rotate(W,phi,c);tau=travel(c,offset);a=tau*(1./c.Q(:));
N=c.nfft;f=(0:N/2)'/(N*c.dt);b=-pi*f+2i*f.*log(max(f,realmin)/c.fref);b(1)=0;
H=exp(b*a');E=fft(eye(c.nt),N,1);G=fromfreq(H,E,c.nt);
op=struct('A',G*Wp,'B',Wp,'dA',{{}},'dB',{{}});
if ~derivatives,return;end
for l=1:4
 da=-tau(:,l)/c.Q(l);dG=fromfreq((b*da').*H,E,c.nt);
 op.dA{l}=dG*Wp;op.dB{l}=zeros(size(W));
end
dW=rotate((-6*z.^2+4*z.^4).*exp(-z.^2),phi,c);
op.dA{5}=G*dW;op.dB{5}=dW;
phase=[0;-1i*ones(N/2-1,1);0;1i*ones(N/2-1,1)];
rot=exp(1i*phi*imag(phase));rot(1)=1;rot(N/2+1)=cos(phi);
drot=phase.*rot;drot(N/2+1)=-sin(phi);
dWp=real(ifft(fft(W,N,1).*drot,[],1));dWp=dWp(1:c.nt,:);
op.dA{6}=G*dWp;op.dB{6}=dWp;
end

function out=rotate(W,phi,c)
if phi==0,out=W;return;end
N=c.nfft;rot=[1;exp(-1i*phi)*ones(N/2-1,1);cos(phi);exp(1i*phi)*ones(N/2-1,1)];
out=real(ifft(fft(W,N,1).*rot,[],1));out=out(1:c.nt,:);
end

function G=fromfreq(H,E,nt)
H(end,:)=real(H(end,:));H=[H;conj(H(end-1:-1:2,:))];
G=real(ifft(H.*E,[],1));G=G(1:nt,:);
end

function tau=travel(c,offset)
one=c;one.Q=ones(1,4);[~,~,~,p]=ray_absorption(c.t,offset,one);
bounds=[0,cumsum(2*c.h./c.vp(1:3)),Inf];t=c.t(:);
tau=max(0,min(t,bounds(2:end))-bounds(1:end-1));tau(~isfinite(tau))=0;
tau=tau./sqrt(max(1-(p.*c.vp).^2,realmin));tau(t<=0,:)=0;tau(t<=0,1)=offset/c.vp(1);
end

function [J,T]=jacobians(op,r,ids,sigma)
J=op.A(:,ids);T=op.B(:,ids);
for l=1:6,J(:,end+1)=op.dA{l}*r;T(:,end+1)=op.dB{l}*r;end
J=J/sigma;
end

function [u,infinite,rnk]=risk(J,T,tol)
% Column scaling is a reparameterisation; G uses the identical scaling.
scale=sqrt(sum(J.^2,1));scale(scale<1e-14)=1;
J=J./scale;T=T./scale;[~,S,V]=svd(J,'econ');v=diag(S);
if isempty(v),u=zeros(size(T,1),1);infinite=u>0;rnk=0;return;end
rnk=sum(v>tol*max(v));Vp=V(:,1:rnk);Z=(T*Vp)./v(1:rnk)';u=sqrt(sum(Z.^2,2));
null=T-(T*Vp)*Vp';infinite=sqrt(sum(null.^2,2))>1e-8*max(sqrt(sum(T.^2,2)),1e-12);
end

function v=wrms(u,flag,mask)
if any(flag(mask)),v=Inf;else,v=sqrt(mean(u(mask).^2));end
end

function [K,diag]=calibration_information(cal,c)
K=zeros(0,6);diag=struct('reduced_residual',NaN,'rank',0);
if isempty(cal),return;end
p=cal.public;p.t=p.t(:);p.vp=reshape(p.vp,1,[]);p.h=reshape(p.h,1,[]);
times=cumsum(2*p.h./p.vp(1:3));mask=false(p.nt,1);
for tt=times,mask=mask|abs(p.t-tt)<=.12;end
resnorm=0;dof=0;h=1e-4;
for j=1:numel(p.offset)
 ni=estimate_trace_snr(cal.D(:,j),p.t);sig=max(ni.sigma,realmin);
 B=p1_layered_calibration_templates(p,p.offset(j),c.Q,c.fpeak);C=B(mask,:)/sig;y=cal.D(mask,j)/sig;
 coef=pinv(C)*y;[U,S,~]=svd(C,'econ');keep=diagvalues(S)>1e-10*norm(C,2);U=U(:,keep);
 J=zeros(sum(mask),6);
 for l=1:5
  qp=c.Q;qm=c.Q;fp=c.fpeak;fm=c.fpeak;
  if l<=4,qp(l)=qp(l)*exp(h);qm(l)=qm(l)*exp(-h);else,fp=fp*exp(h);fm=fm*exp(-h);end
  P=p1_layered_calibration_templates(p,p.offset(j),qp,fp);M=p1_layered_calibration_templates(p,p.offset(j),qm,fm);
  delta=(P(mask,:)-M(mask,:))/(2*h)*coef/sig;J(:,l)=delta-U*(U'*delta);
 end
 % Independent event phase/amplitude nuisance removes a global phase rotation.
 K=[K;J];resnorm=resnorm+norm(y-C*coef)^2;dof=dof+numel(y)-size(U,2); %#ok<AGROW>
end
diag.reduced_residual=sqrt(resnorm/max(dof,1));diag.rank=rank(K,1e-8*max(norm(K,2),1));
end

function v=diagvalues(S)
v=diag(S);
end

function rows=derivative_checks(op,c,off,r,model,condition)
rows={};names={'logQ1','logQ2','logQ3','logQ4','logf','phase_rad','reflectivity','mixed'};
dr=randn(size(r));dr=dr/max(norm(dr),realmin)*max(norm(r),.01);direction=[.13 -.11 .09 -.07 .04 .03];
for l=1:8
 for h=[1e-2 1e-3 1e-4 1e-5]
  delta=zeros(1,6);rr=zeros(size(r));
  if l<=6,delta(l)=1;elseif l==7,rr=dr;else,delta=direction;rr=dr;end
  cp=c;cm=c;cp.Q=c.Q.*exp(h*delta(1:4));cm.Q=c.Q.*exp(-h*delta(1:4));
  cp.fpeak=c.fpeak*exp(h*delta(5));cm.fpeak=c.fpeak*exp(-h*delta(5));
  p=operators(cp,off,h*delta(6),false);m=operators(cm,off,-h*delta(6),false);
  df=(p.A*(r+h*rr)-m.A*(r-h*rr))/(2*h);dt=(p.B*(r+h*rr)-m.B*(r-h*rr))/(2*h);
  jf=op.A*rr;gt=op.B*rr;
  for q=1:6,jf=jf+op.dA{q}*r*delta(q);gt=gt+op.dB{q}*r*delta(q);end
  er=norm(df-jf)/max(norm(jf),1e-12);et=norm(dt-gt)/max(norm(gt),1e-12);
  rows(end+1,:)={model,condition,names{l},h,er,et,norm(dt-gt)}; %#ok<AGROW>
 end
end
end
