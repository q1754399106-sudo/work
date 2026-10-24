function vsp_q_signature()
% POST-HOC exploration E4-1b (decided after E4-1 failed): normalise each level's
% downhole direct-wave spectrum by the same level's surface-hydrophone (source
% signature) spectrum, then repeat spectral-ratio and CFS. Not a registered gate.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(fileparts(root),'bae_q_20261024','source'));
seg=fullfile(root,'raw','ldeo','vsp_segy');[V,dt,H]=read_segy_ibm(fullfile(seg,'308-U1320A_stacks-traces.segy'));
[S,~,HS]=read_segy_ibm(fullfile(seg,'308-U1320A_stacks-surface.segy'));assert(isequal(H.seq,HS.seq));
vq=jsondecode(fileread(fullfile(root,'analysis','vsp_q.json')));keep=find(H.transit_ms>0&H.seq~=14);
nf=1024;f=(0:nf/2)'/(nf*dt);w=taper_pct('taper',100,.2);R=zeros(nf/2+1,numel(keep));srcpk=zeros(1,numel(keep));
for k=1:numel(keep)
 j=keep(k);x=V(:,j);fb=round((H.transit_ms(j)-H.delay_ms(j))/1000/dt)+1;
 s=abs(fft(x(fb-10:fb+89).*w,nf));y=S(:,j);[~,p]=max(abs(y));srcpk(k)=p;
 g=abs(fft(y(max(1,p-30):max(1,p-30)+99).*w,nf));R(:,k)=movmean(s(1:nf/2+1)./max(g(1:nf/2+1),eps),5);
end
depth=arrayfun(@(j)vq.depth_mbsf(find(abs(vq.transit_s-H.transit_ms(j)/1000)<5e-4,1)),keep);tau=H.transit_ms(keep)/1000;
[depth,o]=sort(depth);R=R(:,o);tau=tau(o);band=f>=10.7&f<=119.2;fb_=f(band);
dts=zeros(numel(keep)-1,1);r2=dts;dtau=tau(2:end)'-tau(1);fc=zeros(numel(keep),1);
for k=1:numel(keep)
 p=R(band,k).^2;fc(k)=sum(fb_.*p)/sum(p);if k==1,v1=sum((fb_-fc(1)).^2.*p)/sum(p);continue;end
 y=log(R(band,k)./R(band,1));A=[ones(size(fb_)) fb_];b=A\y;res=y-A*b;dts(k-1)=-b(2)/pi;r2(k-1)=1-sum(res.^2)/sum((y-mean(y)).^2);
end
B=[ones(size(dtau)) dtau]\dts;rng(20261024,'twister');m=numel(dtau);boot=zeros(2000,1);
for b=1:2000,i=randi(m,m,1);c=[ones(m,1) dtau(i)]\dts(i);boot(b)=c(2);end
ci=taper_pct('pct',boot,[5 95]);C=[ones(m,1) dtau]\(fc(1)-fc(2:end));
out=struct('status','POST-HOC exploratory, not a registered gate','Q_SR',1/B(2),'slope_90CI',ci,'median_r2',median(r2), ...
 'Q_CFS',pi*v1/C(2),'centroid_Hz',fc','source_peak_sample',srcpk,'depth_mbsf',depth);
uq_writejson(fullfile(root,'analysis','vsp_q_signature_posthoc.json'),out);disp(jsonencode(out));
end
