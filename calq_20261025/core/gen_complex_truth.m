function model=gen_complex_truth(g)
% Copy of gen_complex_components that also returns the per-trace true attenuation
% a(t) (t* incl. offset path factor), source peak frequency and phase (diagnosis only).
% Original: copy of v1 gen_complex_model that also returns the coherent, multiple and noise
% components (diagnosis only). MATLAB port of complex_fieldlike build_model.py with the
% geometry parameters exposed (g). Quasi-1D column-wise convolution, constant-Q
% attenuation/dispersion, kinematic diffractions and multiples; not full-wave.
rng(g.seed,'twister');nx=99;dx=30;nz=441;dz=5;nt=851;dt=.002;offset=450;
x=(0:nx-1)*dx;z=(0:nz-1)'*dz;t=(0:nt-1)'*dt;[X,Z]=meshgrid(x,z);
fold=g.fold_amp(1)*sin(2*pi*x/2940+g.fold_phase)+g.fold_amp(2)*sin(4*pi*x/2940+.4);
throw1=g.throws(1)*(x>=g.faults(1));throw2=g.throws(2)*(x>=g.faults(2));
strata=Z-fold-(throw1+throw2).*min(max((Z-300)/900,0),1);
horizons=330+[0,cumsum(12+18*rand(1,59))];
layer=1+sum(strata(:)>horizons,2);layer=reshape(layer,size(strata));nlay=numel(horizons)+1;
av=48*randn(nlay,1);ar=.025*randn(nlay,1);aq=.12*randn(nlay,1);av(1)=0;ar(1)=0;aq(1)=0;
base_depth=max(strata-330,0);
vp=g.vp0+.60*base_depth+av(layer).*(1+.15*sin(X/450+(layer-1)*.3));
rho=1.95+.00024*base_depth+ar(layer);
Q=(g.q0+g.q_gradient*base_depth).*exp(aq(layer)+.10*sin(X/700));
for b=1:size(g.steps,1)
 above=strata>g.steps(b,1);vp=vp+g.steps(b,2)*above;rho=rho+.035*above;Q=Q+6*above;
end
channel=((X-g.channel(1))/g.channel(2)).^2+((strata-g.channel(3))/110).^2<1;
vp(channel)=2220+.10*(strata(channel)-g.channel(3));rho(channel)=2.05;Q(channel)=40;
wedge_h=g.wedge_thickness*min(max(1-(x-g.wedge_start)/(g.wedge_tip-g.wedge_start),0),1);
wedge=(strata>=g.wedge_depth)&(strata<g.wedge_depth+wedge_h);
vp=vp+180*wedge;rho=rho+.075*wedge;Q(wedge)=Q(wedge)*.82;
gas=((X-g.gas(1))/330).^2+((strata-g.gas(2))/85).^2<1;
vp=vp-220*gas;rho=rho-.10*gas;Q(gas)=Q(gas)*.52;
het=randn(size(vp));
for it=1:3
 het=(circshift(het,1,1)+2*het+circshift(het,-1,1))/4;het=(circshift(het,1,2)+2*het+circshift(het,-1,2))/4;
end
mask=strata>330;vp=vp+12*het.*mask;rho=rho+.006*het.*mask;
vp=min(max(vp,1750),3650);rho=min(max(rho,1.85),2.65);Q=min(max(Q,18),160);
vs=vp.*(.46+.015*sin(strata/370)).*(1+.04*gas);vs(strata<=330)=874;
N=4096;f=(0:N/2)'/(N*dt);bq=-pi*f+2i*f.*log(max(f,realmin)/100);bq(1)=0;Ef=exp(-2i*pi*f*t');
att=@(s,a)attenuate(s,a,bq,Ef,N,nt);
primary=zeros(nt,nx);atten=primary;multiple=primary;cumtime=zeros(nz,nx);atrue=primary;fptrue=zeros(1,nx);phtrue=fptrue;
for j=1:nx
 xx=x(j);travel=[0;cumsum(2*dz./((vp(1:end-1,j)+vp(2:end,j))/2))];cumtime(:,j)=travel;
 tau=diff(travel);vr=.5*(vp(1:end-1,j)+vp(2:end,j));rr=.5*(rho(1:end-1,j)+rho(2:end,j));vsr=.5*(vs(1:end-1,j)+vs(2:end,j));
 A=.5*(diff(vp(:,j))./vr+diff(rho(:,j))./rr);
 B=.5*diff(vp(:,j))./vr-2*(vsr./vr).^2.*(diff(rho(:,j))./rr+2*diff(vs(:,j))./vsr);
 events=travel(2:end);vrms=sqrt(cumsum(vr.^2.*tau)./events);
 sin2=offset^2./(offset^2+(vrms.*events).^2);coeff=A+B.*sin2;
 fpeak=36*(1+.07*sin(2*pi*xx/2940+.3));phase=24+7*sin(xx/560);
 zz=pi*fpeak*(t-events');W=(1-2*zz.^2).*exp(-zz.^2);clean=W*coeff;
 for d=1:size(g.diffractors,1)
  c=g.diffractors(d,:);tt=sqrt(c(2)^2+4*(xx-c(1))^2/2400^2);u=pi*fpeak*(t-tt);
  clean=clean+c(3)/(1+abs(xx-c(1))/300)*(1-2*u.^2).*exp(-u.^2);
 end
 clean=phase_source(clean,phase,t,N,nt);
 aqcum=[0;cumsum(tau./sqrt(Q(1:end-1,j).*Q(2:end,j)))];a0=interp1(travel,aqcum,min(t,travel(end)),'linear');
 vtime=interp1(events,vrms,t,'linear',NaN);vtime(t<events(1))=vp(1,j);vtime(t>events(end))=vrms(end);
 a=a0.*sqrt(1+(offset./max(vtime.*t,vp(1,j)*.02)).^2);
 primary(:,j)=clean;atten(:,j)=att(clean,a);atrue(:,j)=a;fptrue(j)=fpeak;phtrue(j)=phase;
 dirty=.14*delay(clean,.18+.025*sin(xx/600),t)+.08*delay(clean,.36,t);multiple(:,j)=att(dirty,a);
end
% Independent coarse four-layer calibration, same source family.
edges=[0 330 780 1260 2205];cvp=zeros(1,4);cq=cvp;
for L=1:4,k=z>=edges(L)&z<edges(L+1);cvp(L)=mean(vp(k,:),'all');cq(L)=1/mean(1./Q(k,:),'all');end
ch=diff(edges(1:4));ctime=cumsum(2*ch./cvp(1:3));coff=(0:10)*150;
public=struct('dt',dt,'nt',nt,'t',t,'offset',coff,'vp',cvp,'h',ch,'fref',100,'nfft',N,'lambda_x',.001, ...
 'lambda_fixed',.001,'lambda_grid',logspace(-7,0,41),'cov_floor',1e-6,'nominal_frequency',36);
cal=zeros(nt,11);calclean=cal;rng(g.calibration_seed,'twister');
cfg=struct('dt',dt,'nt',nt,'h',ch,'vp',cvp,'Q',cq);
for j=1:11
 [a0,ax,~,p]=ray_absorption(t,coff(j),cfg);a=a0+ax;a(t<=0)=coff(j)/(cvp(1)*cq(1));
 pp=interp1(t,p,ctime);am=[.10,-.045,.075]+[-.045,.025,-.030].*(pp.*cvp(1:3)).^2;
 u=pi*36*(t-ctime);s=((1-2*u.^2).*exp(-u.^2))*am';s=phase_source(s,24+7*sin(1470/560),t,N,nt);s=att(s,a);
 noise=randn(nt,1);cal(:,j)=s+noise/sqrt(mean(noise.^2))*sqrt(mean(s.^2))*10^(-18/20);calclean(:,j)=s;
end
live=true(1,nx);live(g.dead)=false;datasets=struct('name',{},'D',{},'target',{},'coherent',{},'multiple',{},'noise',{},'primary_att',{});
for o=1:numel(g.observations)
 cfgo=g.observations(o);rng(cfgo.seed,'twister');gain=exp(.09*randn(1,nx));stat=cfgo.static_std_ms*randn(1,nx);
 target=zeros(nt,nx);base=target;arms=zeros(1,nx);
 for j=1:nx
  target(:,j)=gain(j)*delay(primary(:,j),stat(j)*1e-3,t);
  base(:,j)=gain(j)*delay(atten(:,j)+cfgo.multiple_scale*multiple(:,j),stat(j)*1e-3,t);
  arms(j)=sqrt(mean((gain(j)*delay(atten(:,j),stat(j)*1e-3,t)).^2));
 end
 noise=randn(nt,nx);for k=2:nt,noise(k,:)=.55*noise(k-1,:)+sqrt(1-.55^2)*noise(k,:);end
 noise=(circshift(noise,1,2)+2*noise+circshift(noise,-1,2))/4;
 noise=noise./sqrt(mean(noise.^2,1)).*arms*10^(-cfgo.snr_db/20);
 coherent=zeros(nt,nx);
 for j=1:nx
  u=t-(.18+abs(x(j)-g.groundroll_x)/1850);ground=sin(2*pi*(10*u+8*u.^2)).*exp(-(u/.095).^2);
  v=t-(g.linear_t0+g.linear_slope*x(j));lin=cos(2*pi*29*v).*exp(-(v/.045).^2);coherent(:,j)=ground+.35*lin;
 end
 coherent=coherent*cfgo.coherent_fraction*norm(base,'fro')/norm(coherent,'fro');bursts=zeros(nt,nx);
 for j=g.burst_traces,idx=randi([171 nt-50],4,1);bursts(idx,j)=sign(randn(4,1))*6*arms(j);end
 D=base+noise+coherent+bursts;D(:,~live)=0;
 mult=zeros(nt,nx);patt=mult;for j=1:nx,mult(:,j)=gain(j)*delay(cfgo.multiple_scale*multiple(:,j),stat(j)*1e-3,t);patt(:,j)=gain(j)*delay(atten(:,j),stat(j)*1e-3,t);end
 datasets(o)=struct('name',cfgo.name,'D',D,'target',target,'coherent',coherent,'multiple',mult,'noise',noise+bursts,'primary_att',patt);
end
model=struct('x',x,'t',t,'z',z,'public',public,'calibrationD',cal,'calibration_clean',calclean,'datasets',datasets,'live',live, ...
 'truth',struct('vp',vp,'Q',Q,'strata',strata,'cumtime',cumtime,'wedge',wedge,'channel',channel,'gas',gas,'calibration_Q',cq,'a',atrue,'fpeak',fptrue,'phase',phtrue),'geometry',g);
end

function y=attenuate(s,a,b,E,N,nt)
Y=sum(exp(b*a').*E.*s',2);Y(end)=real(Y(end));
y=real(ifft([Y;conj(Y(end-1:-1:2))]));y=y(1:nt);
end

function w=phase_source(s,phase,t,N,nt)
F=fft(s,N);F=F(1:N/2+1);rot=exp(-1i*deg2rad(phase))*ones(size(F));rot(1)=1;rot(end)=cos(deg2rad(phase));
G=F.*rot;w=real(ifft([G;conj(G(end-1:-1:2))]));w=w(1:nt);
w=w-.16*delay(w,.012,t)+.04*delay(w,.025,t);
end

function y=delay(s,shift,t)
y=interp1(t,s,t-shift,'linear',0);
end
