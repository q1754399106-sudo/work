function vsp_q()
% E6-1 (PLAN_ZH.md): interval Q of segments RS and NI from level-to-level spectral ratios (reference =
% shallowest level of the segment) with CFS and G1; Q_top from the ghost-corrected near-field hydrophone of
% the RS reference level with G1b; the 4-layer model Q_VSP. Spectral computations follow experiment 5.
P=setup_paths();root=P.root;raw=fullfile(root,'raw','vsp_velocity');
geo=jsondecode(fileread(fullfile(root,'probe','vsp_geometry.json')));
rs=asrow(geo.RS);rs=rs([rs.md]<=P.rs_md_max);ni=asrow(geo.NI);
[RS,f,w,nf]=segment(P,fullfile(raw,'VSPRS_COMPUTED_1.SEGY'),rs,'RS');
NI=segment(P,fullfile(raw,'VSPNI_COMPUTED_1.SEGY'),ni,'NI');
% Source-referenced t* to the RS reference level (near-field hydrophone, ghost-corrected far field).
[N,dtn,Hn]=read_segy(fullfile(raw,'VSPRS_RAW_1.SEGY'));rec=hword(Hn,9,'int32');
recs=RS.ref_level.nfh_records(:)';cols=arrayfun(@(r)find(rec==r,1),recs);assert(numel(cols)==numel(recs));
dt=RS.dt;step=round(dt/dtn);segs=zeros(100,numel(cols));
for q=1:numel(cols)
 x=decimate_fft(N(:,cols(q)),step);[~,pk]=max(abs(x));segs(:,q)=x(pk-10:pk+89);
end
ghost=abs(1-exp(-1i*2*pi*f*2*5.63/1480));sref=RS.sig(:,RS.ref);good=RS.good;band=f(good);
[bt,r2top]=top_fit(segs,band,good,sref,w,nf,ghost);tstar_top=-bt(2)/pi;
rng(P.seed+1,'twister');bb=zeros(2000,1);nr=numel(cols);
for b=1:2000,i=randi(nr,nr,1);c=top_fit(segs(:,i),band,good,sref,w,nf,ghost);bb(b)=-c(2)/pi;end
ci_top=taper_pct('pct',bb,[5 95]);G1b=ci_top(1)>0&&r2top>=.5;
Qtop_est=RS.fb(RS.ref)/tstar_top;
% Segment fallbacks (PLAN_ZH.md "Q 模型").
okRS=isfinite(RS.Q)&&RS.Q>0;okNI=isfinite(NI.Q)&&NI.Q>0;Qrs=RS.Q;Qni=NI.Q;notes={};
if ~okRS&&okNI,Qrs=Qni;notes{end+1}='Q_RS invalid, replaced by Q_NI';end
if ~okNI&&okRS,Qni=Qrs;notes{end+1}='Q_NI invalid, replaced by Q_RS';end
executable=okRS||okNI;
Qtop=Qtop_est;if ~G1b,Qtop=Qrs;notes{end+1}='G1b failed: Q_top = Q_RS';end
% Boundaries (vertical two-way time) and velocities from the contractor's combined time-depth table.
cs=read_checkshot(fullfile(raw,'TZV_TIME_CHECKSHOT_3.ASC'));
mdb=[RS.ref_level.md NI.ref_level.md max([ni.md])];
tb=interp1(cs(:,1),cs(:,3),mdb)/1000;zb=interp1(cs(:,1),cs(:,2),mdb);
vp=[zb(1)/(tb(1)/2) diff(zb)./(diff(tb)/2)];vp(4)=vp(3);h=vp(1:3).*diff([0 tb])/2;
models=struct('VSP',[Qtop Qrs Qni Qni]);
G1=struct('RS',RS.gate,'NI',NI.gate);G1.all_pass=RS.gate.G1_pass&&NI.gate.G1_pass;
out=struct('RS',rmfield(RS,{'sig','noi','good','ref_level','dt'}),'NI',rmfield(NI,{'sig','noi','good','ref_level','dt'}), ...
 'G1',G1,'tstar_top',tstar_top,'tstar_top_90CI',ci_top,'top_r2',r2top,'Q_top_estimate',Qtop_est,'G1b_pass',G1b, ...
 'Q_top_used',Qtop,'nfh_records',recs,'ghost_depth_m',5.63,'bounds_twt',tb,'bounds_md',mdb,'vp',vp,'h',h, ...
 'Q_models',models,'V_executable',executable,'notes',{notes});
if ~isfolder(fullfile(root,'analysis')),mkdir(fullfile(root,'analysis'));end
uq_writejson(fullfile(root,'analysis','vsp_q.json'),out);
disp(jsonencode(struct('G1',G1,'Q_RS',RS.Q,'Q_RS_90CI',RS.Q_90CI,'Q_NI',NI.Q,'Q_NI_90CI',NI.Q_90CI, ...
 'Q_top_estimate',Qtop_est,'G1b_pass',G1b,'Q_VSP',models.VSP,'bounds_twt',tb,'notes',{notes})));
sigRS=RS.sig;noiRS=RS.noi;sigNI=NI.sig;noiNI=NI.noi;
save(fullfile(root,'analysis','vsp_spectra.mat'),'f','sigRS','noiRS','sigNI','noiNI','segs','ghost');
end

function [S,f,w,nf]=segment(P,file,lv,name)
% Edit, spectra, band, spectral-ratio Q with bootstrap CI, CFS Q (power-spectrum moments, 2*pi) and G1.
[~,o]=sort([lv.md]);lv=lv(o);fb=[lv.fb_ms]'/1000;
[V,dt]=read_segy(file,[lv.trace]);
dead=sqrt(mean(V.^2,1))'==0;steps=diff(fb);med=median(steps);keep=~dead;last=find(keep,1);
for k=last+1:numel(fb)
 if ~keep(k),continue;end
 d=fb(k)-fb(last);if d<0||d>3*med,keep(k)=false;else,last=k;end
end
dropped=[lv(~keep).trace];lv=lv(keep);V=V(:,keep);fb=fb(keep);nl=numel(lv);
nf=1024;f=(0:nf/2)'/(nf*dt);w=taper_pct('taper',100,.2);sig=zeros(numel(f),nl);noi=sig;
for k=1:nl
 i0=round(fb(k)/dt)+1;sig(:,k)=spec(V(i0-10:i0+89,k),w,nf);noi(:,k)=spec(V(i0-210:i0-111,k),w,nf);
end
good=all(sig>=3*noi,2)&f>=10&f<=120;band=f(good);ok_band=numel(band)>=2&&(max(band)-min(band))>=20;
assert(numel(band)>=2,'%s: no usable frequency band',name);
ref=1;oth=2:nl;dts=zeros(numel(oth),1);r2=dts;dtau=fb(oth)-fb(ref);fc=zeros(nl,1);
for k=1:nl,p=sig(good,k).^2;fc(k)=sum(band.*p)/sum(p);end
pr=sig(good,ref).^2;var_ref=sum((band-fc(ref)).^2.*pr)/sum(pr);
for i=1:numel(oth),[b,r2(i)]=linfit(band,log(sig(good,oth(i))./sig(good,ref)));dts(i)=-b(2)/pi;end
[B,reg_r2]=linfit(dtau,dts);slope=B(2);Q=1/slope;
rng(P.seed,'twister');boot=zeros(2000,1);m=numel(dtau);
for b=1:2000,i=randi(m,m,1);if numel(unique(dtau(i)))<2,boot(b)=NaN;continue;end;c=[ones(m,1) dtau(i)]\dts(i);boot(b)=c(2);end
boot=boot(isfinite(boot));ci=taper_pct('pct',boot,[5 95]);Qci=sort(1./ci);if ci(1)<=0,Qci=[1/ci(2) Inf];end
C=[ones(m,1) dtau]\(fc(ref)-fc(oth));Qcfs=2*pi*var_ref/C(2);
g=struct('median_r2',median(r2),'median_r2_ok',median(r2)>=.5,'slope_ci_positive',ci(1)>0, ...
 'cfs_ratio',Qcfs/Q,'cfs_ok',Qcfs/Q>=.5&&Qcfs/Q<=2,'band_ok',ok_band);
g.G1_pass=g.median_r2_ok&&g.slope_ci_positive&&g.cfs_ok&&g.band_ok;
S=struct('levels_used',nl,'dropped_traces',dropped,'reference_trace',lv(ref).trace,'ref',ref,'ref_level',lv(ref), ...
 'band_Hz',[min(band) max(band)],'band_points',numel(band),'delta_tau',dtau','delta_tstar',dts','level_r2',r2', ...
 'regression_r2',reg_r2,'Q',Q,'Q_90CI',Qci,'slope_90CI',ci,'Q_CFS',Qcfs,'centroid_Hz',fc','gate',g, ...
 'fb',fb,'dt',dt,'sig',sig,'noi',noi,'good',good);
end

function A=spec(x,w,nf)
A=abs(fft(x(:).*w,nf));A=movmean(A(1:nf/2+1),5);
end
function [b,r2]=linfit(x,y)
A=[ones(size(x)) x];b=A\y;res=y-A*b;r2=1-sum(res.^2)/sum((y-mean(y)).^2);
end
function [b,r2]=top_fit(segs,band,good,sref,w,nf,ghost)
S=spec(median(segs,2),w,nf).*ghost;[b,r2]=linfit(band,log(sref(good)./S(good)));
end
function y=decimate_fft(x,step)
% Anti-alias: cosine roll-off from 0.8 to 1.0 of the new Nyquist, then keep every step-th sample.
n=numel(x);X=fft(x);f=[0:floor(n/2) -ceil(n/2)+1:-1]'/n;fn=.5/step;a=abs(f)/fn;
H=ones(n,1);r=a>.8&a<1;H(r)=.5*(1+cos(pi*(a(r)-.8)/.2));H(a>=1)=0;y=real(ifft(X.*H));y=y(1:step:end);
end
function c=read_checkshot(path)
% Rows of three numbers (MD RKB m, TVD MSL m, TWT ms) from a contractor time-depth listing.
L=splitlines(fileread(path));c=zeros(0,3);
for k=1:numel(L),v=sscanf(L{k},'%f')';if numel(v)==3,c(end+1,:)=v;end,end %#ok<AGROW>
[~,u]=unique(c(:,1));c=c(u,:);
end
function s=asrow(s),if iscell(s),s=[s{:}];end,s=s(:)';end
