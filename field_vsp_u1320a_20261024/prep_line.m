function prep_line(line,well_trace)
% Data-driven preparation of a 99-trace window around the well (PLAN_ZH.md).
% Needs analysis/vsp_q.json (checkshot and Q_SR). Writes inputs/<line>.mat.
root=fileparts(mfilename('fullpath'));bae=fullfile(fileparts(root),'bae_q_20261024');
addpath(fullfile(root,'core'),fullfile(bae,'core'),fullfile(bae,'source'));
vq=jsondecode(fileread(fullfile(root,'analysis','vsp_q.json')));
[X,dt,H]=read_segy_ibm(fullfile(root,'raw','namss','segy',[line '.segy']),well_trace-49,99);
t=(0:size(X,1)-1)'*dt;win=t>=1.9-1e-9&t<=3.7+1e-9;tw=t(win);
% Seafloor: first strong sample in 1.8-2.2 s, refined to envelope peak within +-20 ms.
env=abs(X+1i*taper_pct('hilbert',X));sfi=zeros(1,99);
for j=1:99
 r=find(t>=1.8&t<=2.2);k=r(find(abs(X(r,j))>.3*max(abs(X(r,j))),1));
 rr=max(1,k-5):k+5;[~,m]=max(env(rr,j));sfi(j)=rr(m);
end
Tsf=median(t(sfi));
% Seafloor wavelet: aligned average, +-48 ms.
L=12;Wv=zeros(2*L+1,1);for j=1:99,Wv=Wv+X(sfi(j)-L:sfi(j)+L,j)/99;end
nf=512;f=(0:nf/2)'/(nf*dt);A=abs(fft(Wv.*taper_pct('taper',2*L+1,.3),nf));A=A(1:nf/2+1);band=f>=5&f<=100;
grid=10:.25:80;mis=zeros(size(grid));
for q=1:numel(grid),r=(f(band).^2/grid(q)^2).*exp(-f(band).^2/grid(q)^2);a=r\A(band);mis(q)=norm(A(band)-a*r);end
[~,q]=min(mis);fp=grid(q);
tl=(-L:L)'*dt;z=pi*fp*tl;ric=(1-2*z.^2).*exp(-z.^2);Hw=taper_pct('hilbert',Wv);best=-Inf;phi=0;
for ph=-179:180
 y=cosd(ph)*Wv-sind(ph)*Hw;c=(y'*ric)/(norm(y)*norm(ric));if c>best,best=c;phi=ph;end
end
Xr=cosd(phi)*X-sind(phi)*taper_pct('hilbert',X);D=Xr(win,:);
% Velocity model from VSP checkshot (one-way transit, mbsf).
dep=vq.depth_mbsf(:);tau=vq.transit_s(:);s2=dep<=150;s3=dep>=147;
p2=[ones(nnz(s2),1) tau(s2)]\dep(s2);p3=[ones(nnz(s3),1) tau(s3)]\dep(s3);vp2=p2(2);vp3=p3(2);
cm=struct('dt',dt,'nt',numel(tw),'nfft',2048,'fref',100,'fpeak',fp,'t',tw,'vp',[1500 vp2 vp3 1.1*vp3], ...
 'h',[1500*Tsf/2 150 150],'Q',[1e6 vq.Q_SR vq.Q_SR vq.Q_SR],'solver_offset',0,'lambda_x',.001,'lambda_fixed',.001);
sigma=zeros(1,99);for j=1:99,ni=estimate_trace_snr(D(:,j),tw,struct('method','haar'));sigma(j)=ni.sigma;end
meta=struct('line',line,'well_trace',well_trace,'traces',[well_trace-49 well_trace+49],'sp',H.sp([1 50 99]), ...
 'seafloor_twt',Tsf,'seafloor_twt_range',[min(t(sfi)) max(t(sfi))],'fpeak',fp,'phase_deg',phi,'phase_corr',best, ...
 'vp',cm.vp,'h',cm.h,'Q_SR',vq.Q_SR,'sigma_median',median(sigma),'note','Q of layers 2-4 = VSP Q_SR; layer-4 velocity and Q extrapolated');
if ~isfolder(fullfile(root,'inputs')),mkdir(fullfile(root,'inputs'));end
save(fullfile(root,'inputs',[line '.mat']),'D','cm','sigma','meta','Wv','-v7');uq_writejson(fullfile(root,'inputs',[line '_meta.json']),meta);
disp(jsonencode(meta));
end
