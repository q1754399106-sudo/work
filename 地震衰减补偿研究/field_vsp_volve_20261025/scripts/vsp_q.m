function vsp_q()
% E5-1 (PLAN_ZH.md): Q_int from level-to-level spectral ratios (reference level 191) with CFS and G1;
% Q_top from the ghost-corrected near-field hydrophone reference with G1b; the two 4-layer Q models.
P=setup_paths();root=P.root;nl=P.levels;
lv=jsondecode(fileread(fullfile(root,'probe','vsp_geometry.json'))).levels;if iscell(lv),lv=[lv{:}];end
lv=lv(1:nl);fb=[lv.first_break_s]';tvd=[lv.tvd_msl]';
[V,dt]=read_segy(fullfile(root,'raw','vsp_computed','VSPNI_COMPUTED_1.SEGY'),1:nl);
nf=1024;f=(0:nf/2)'/(nf*dt);w=taper_pct('taper',100,.2);sig=zeros(numel(f),nl);noi=sig;
for k=1:nl
 i0=round(fb(k)/dt)+1;sig(:,k)=spec(V(i0-10:i0+89,k),w,nf);noi(:,k)=spec(V(i0-210:i0-111,k),w,nf);
end
good=all(sig>=3*noi,2)&f>=10&f<=120;band=f(good);ok_band=numel(band)>=2&&(max(band)-min(band))>=20;
assert(numel(band)>=2,'no usable frequency band');
% Level-to-level spectral ratios against level 191 (shallowest usable).
ref=nl;oth=1:nl-1;dts=zeros(numel(oth),1);r2=dts;dtau=fb(oth)-fb(ref);fc=zeros(nl,1);
for k=1:nl
 p=sig(good,k).^2;fc(k)=sum(band.*p)/sum(p);
end
pr=sig(good,ref).^2;var_ref=sum((band-fc(ref)).^2.*pr)/sum(pr);
for i=1:numel(oth)
 [b,r2(i)]=linfit(band,log(sig(good,oth(i))./sig(good,ref)));dts(i)=-b(2)/pi;
end
B=[ones(size(dtau)) dtau]\dts;slope=B(2);Qint=1/slope;
rng(P.seed,'twister');boot=zeros(2000,1);m=numel(dtau);
for b=1:2000,i=randi(m,m,1);if numel(unique(dtau(i)))<2,boot(b)=NaN;continue;end;c=[ones(m,1) dtau(i)]\dts(i);boot(b)=c(2);end
boot=boot(isfinite(boot));ci_slope=taper_pct('pct',boot,[5 95]);Qci=sort(1./ci_slope);if ci_slope(1)<=0,Qci=[1/ci_slope(2) Inf];end
C=[ones(m,1) dtau]\(fc(ref)-fc(oth));Qcfs=pi*var_ref/C(2);
g=struct('median_r2',median(r2),'median_r2_ok',median(r2)>=.5,'slope_ci_positive',ci_slope(1)>0, ...
 'cfs_int_ratio',Qcfs/Qint,'cfs_int_ok',Qcfs/Qint>=.5&&Qcfs/Qint<=2,'band_ok',ok_band);
g.G1_pass=g.median_r2_ok&&g.slope_ci_positive&&g.cfs_int_ok&&g.band_ok;
% Source-referenced t* to level 191 from the near-field hydrophone (ghost-corrected far-field estimate).
[N,dtn,Hn]=read_segy(fullfile(root,'raw','vsp_raw','VSPNI_RAW_1.SEGY'));rec=hword(Hn,9,'int32');
recs=lv(ref).nfh_records(:)';cols=arrayfun(@(r)find(rec==r,1),recs);assert(numel(cols)==numel(recs));
step=round(dt/dtn);segs=zeros(100,numel(cols));
for q=1:numel(cols)
 x=decimate_fft(N(:,cols(q)),step);[~,pk]=max(abs(x));segs(:,q)=x(pk-10:pk+89);
end
ghost=abs(1-exp(-1i*2*pi*f*2*5.0/1480));
[bt,r2top]=top_fit(segs,band,good,sig(:,ref),w,nf,ghost);tstar_top=-bt(2)/pi;
bb=zeros(2000,1);nr=numel(cols);
for b=1:2000,i=randi(nr,nr,1);c=top_fit(segs(:,i),band,good,sig(:,ref),w,nf,ghost);bb(b)=-c(2)/pi;end
ci_top=taper_pct('pct',bb,[5 95]);G1b=ci_top(1)>0&&r2top>=.5;
Qtop_est=fb(ref)/tstar_top;Qtop=Qtop_est;if ~G1b,Qtop=Qint;end
% Four-layer models: two-way boundaries [T191 2.415 2.595]; h from VSP interval velocities.
bounds=[2*fb(ref) 2.415 2.595];vp=zeros(1,4);edges=[0 bounds Inf];
vp(1)=tvd(ref)/fb(ref);
for L=2:4
 s=2*fb>=edges(L)&2*fb<=min(edges(L+1),2*max(fb));p=[ones(nnz(s),1) fb(s)]\tvd(s);vp(L)=p(2);
end
h=vp(1:3).*diff([0 bounds])/2;
models=struct('VSP',[Qtop Qint Qint Qint],'READ',[46.22496 300 49.81825 300]);
out=struct('levels_used',nl,'reference_level',ref,'band_Hz',[min(band) max(band)],'band_points',numel(band), ...
 'delta_tau',dtau','delta_tstar',dts','level_r2',r2','Q_int',Qint,'Q_int_90CI',Qci,'slope_90CI',ci_slope, ...
 'Q_CFS',Qcfs,'centroid_Hz',fc','gate',g,'tstar_top',tstar_top,'tstar_top_90CI',ci_top,'top_r2',r2top, ...
 'Q_top_estimate',Qtop_est,'G1b_pass',G1b,'Q_top_used',Qtop,'nfh_records',recs,'bounds_twt',bounds,'vp',vp,'h',h, ...
 'Q_models',models);
if ~isfolder(fullfile(root,'analysis')),mkdir(fullfile(root,'analysis'));end
uq_writejson(fullfile(root,'analysis','vsp_q.json'),out);disp(jsonencode(rmfield(out,{'delta_tau','delta_tstar','level_r2','centroid_Hz'})));
save(fullfile(root,'analysis','vsp_spectra.mat'),'f','sig','noi','good','fb','tvd','segs','ghost');
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
