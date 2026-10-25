function vsp_q()
% E4-1: spectral-ratio and centroid-frequency-shift Q from U1320A VSP stacks (PLAN_ZH.md).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(fileparts(root),'bae_q_20261024','source'));
[V,dt,H]=read_segy_ibm(fullfile(root,'raw','ldeo','vsp_segy','308-U1320A_stacks-traces.segy'));
lv=vsp_table(root);seq=H.seq;keep=find(H.transit_ms>0&seq~=14);
depth=arrayfun(@(s)lv.mbsf(lv.seq==s&lv.transit>0),seq(keep));tau=H.transit_ms(keep)/1000;
[depth,o]=sort(depth);keep=keep(o);tau=tau(o);
nf=1024;f=(0:nf/2)'/(nf*dt);sig=zeros(numel(f),numel(keep));noi=sig;
for k=1:numel(keep)
 x=V(:,keep(k));fb=round((H.transit_ms(keep(k))-H.delay_ms(keep(k)))/1000/dt)+1;
 s=x(fb-10:fb+89);n=x(fb-210:fb-111);w=taper_pct('taper',100,.2);
 S=abs(fft(s.*w,nf));N=abs(fft(n.*w,nf));sig(:,k)=movmean(S(1:nf/2+1),5);noi(:,k)=movmean(N(1:nf/2+1),5);
end
good=all(sig>=3*noi,2)&f>=10&f<=120;band=f(good);
ok_band=numel(band)>=2&&(max(band)-min(band))>=20;
ref=1;dts=zeros(numel(keep)-1,1);r2=dts;dtau=tau(2:end)'-tau(ref);fc=zeros(numel(keep),1);
for k=1:numel(keep)
 p=sig(good,k).^2;fc(k)=sum(band.*p)/sum(p);
 if k==ref,var_ref=sum((band-fc(k)).^2.*p)/sum(p);continue;end
 y=log(sig(good,k)./sig(good,ref));A=[ones(size(band)) band];b=A\y;res=y-A*b;
 dts(k-1)=-b(2)/pi;r2(k-1)=1-sum(res.^2)/sum((y-mean(y)).^2);
end
B=[ones(size(dtau)) dtau]\dts;slope=B(2);Qsr=1/slope;
rng(20261024,'twister');boot=zeros(2000,1);m=numel(dtau);
for b=1:2000,i=randi(m,m,1);if numel(unique(dtau(i)))<2,boot(b)=NaN;continue;end;c=[ones(m,1) dtau(i)]\dts(i);boot(b)=c(2);end
boot=boot(isfinite(boot));ci_slope=taper_pct('pct',boot,[5 95]);Qci=sort(1./ci_slope);if ci_slope(1)<=0,Qci(2)=Inf;end
C=[ones(m,1) dtau]\(fc(ref)-fc(2:end));Qcfs=pi*var_ref/C(2);
g=struct('median_r2',median(r2),'median_r2_ok',median(r2)>=.5,'slope_ci_positive',ci_slope(1)>0, ...
 'cfs_sr_ratio',Qcfs/Qsr,'cfs_sr_ok',Qcfs/Qsr>=.5&&Qcfs/Qsr<=2,'band_ok',ok_band);
g.G1_pass=g.median_r2_ok&&g.slope_ci_positive&&g.cfs_sr_ok&&g.band_ok;
out=struct('levels_used',numel(keep),'depth_mbsf',depth,'transit_s',tau,'band_Hz',[min(band) max(band)], ...
 'delta_tstar',dts','delta_tau',dtau','level_r2',r2','Q_SR',Qsr,'Q_SR_90CI',Qci,'slope_90CI',ci_slope, ...
 'Q_CFS',Qcfs,'centroid_Hz',fc','gate',g);
if ~isfolder(fullfile(root,'analysis')),mkdir(fullfile(root,'analysis'));end
uq_writejson(fullfile(root,'analysis','vsp_q.json'),out);disp(jsonencode(out));
save(fullfile(root,'analysis','vsp_spectra.mat'),'f','sig','noi','depth','tau','good');
end

function lv=vsp_table(root)
txt=splitlines(fileread(fullfile(root,'raw','ldeo','vsp_ascii','308-U1320A_stacks-summary.txt')));
seq=[];mbsf=[];tr=[];
for i=1:numel(txt)
 p=strsplit(strtrim(txt{i}));
 if numel(p)==6&&~isnan(str2double(p{1}))&&~isnan(str2double(p{2}))
  seq(end+1)=str2double(p{1});mbsf(end+1)=str2double(p{4})-1478;tr(end+1)=str2double(p{6}); %#ok<AGROW>
 end
end
lv=struct('seq',seq,'mbsf',mbsf,'transit',tr);
end
