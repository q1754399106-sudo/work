function diag_source()
% Is the calibration bias caused by source mismatch (ghost / phase) in the
% calibration templates? Rebuild the noise-free calibration gather with the
% generator's code, switching ghost and phase, and run the frozen calibration.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));out=struct();
for v={'V3','V5'}
 TT=load(fullfile(root,'dev',[v{1} '_truth.mat']),'g');model=gen_complex_truth(TT.g);p=model.public;cq=model.truth.calibration_Q;
 ref=calib_gather(p,cq,24+7*sin(1470/560),true);
 fprintf('%s rebuild vs generator noise-free calibration: %.3g\n',v{1},norm(ref-model.calibration_clean,'fro')/norm(model.calibration_clean,'fro'));
 variants={'ghost+phase',24+7*sin(1470/560),true;'phase only',24+7*sin(1470/560),false;'ghost only',0,true;'neither',0,false};
 for k=1:size(variants,1)
  D=calib_gather(p,cq,variants{k,2},variants{k,3});e=p1_prepare_candidate(D,p.offset,p,policy);
  fprintf('%s %-12s Q %s  ln(Q/Qtrue) %s\n',v{1},variants{k,1},mat2str(e.layer_Q,4),mat2str(log(e.layer_Q./cq),3));
  out.(v{1}).(matlab.lang.makeValidName(variants{k,1}))=struct('Q',e.layer_Q,'log_bias',log(e.layer_Q./cq),'true',cq);
 end
end
uq_writejson(fullfile(root,'analysis','diag_source.json'),out);
end

function cal=calib_gather(public,cq,phase,ghost)
% Calibration section of gen_complex_truth.m (noise-free), with switchable phase and ghost.
t=public.t;nt=public.nt;N=public.nfft;dt=public.dt;cvp=public.vp;ch=public.h;coff=public.offset;
f=(0:N/2)'/(N*dt);bq=-pi*f+2i*f.*log(max(f,realmin)/100);bq(1)=0;Ef=exp(-2i*pi*f*t');
cfg=struct('dt',dt,'nt',nt,'h',ch,'vp',cvp,'Q',cq);ctime=cumsum(2*ch./cvp(1:3));cal=zeros(nt,11);
for j=1:11
 [a0,ax,~,pr]=ray_absorption(t,coff(j),cfg);a=a0+ax;a(t<=0)=coff(j)/(cvp(1)*cq(1));
 pp=interp1(t,pr,ctime);am=[.10,-.045,.075]+[-.045,.025,-.030].*(pp.*cvp(1:3)).^2;
 u=pi*36*(t-ctime);s=((1-2*u.^2).*exp(-u.^2))*am';
 F=fft(s,N);F=F(1:N/2+1);rot=exp(-1i*deg2rad(phase))*ones(size(F));rot(1)=1;rot(end)=cos(deg2rad(phase));
 Gf=F.*rot;w=real(ifft([Gf;conj(Gf(end-1:-1:2))]));w=w(1:nt);
 if ghost,w=w-.16*interp1(t,w,t-.012,'linear',0)+.04*interp1(t,w,t-.025,'linear',0);end
 Y=sum(exp(bq*a').*Ef.*w',2);Y(end)=real(Y(end));y=real(ifft([Y;conj(Y(end-1:-1:2))]));cal(:,j)=y(1:nt);
end
end
