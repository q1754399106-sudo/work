function check_octave(j)
% Octave round trip on archived data (no new data): geology H7 of fd2d_confirm_20261109 (MATLAB), shots 0 and j in Octave
% with octave_compat_v2, processed with process_zo, against the archived clean trace j of models/H7.mat.
here=fileparts(mfilename('fullpath'));addpath(fileparts(here));
P=fd_setup_h712();sp=P.models(1);t0=tic;Gm=fd_geology(sp,P);
A=load(fullfile(P.root,'..','fd2d_confirm_20261109','models','H7.mat'));
fprintf('geology: %d fine layers, upscaled Q %s (archived %s), max |dQ| %.2e  (%.0fs)\n',numel(Gm.base.v),mat2str(Gm.config.Q,6),mat2str(A.config.Q,6),max(abs(Gm.config.Q(:)-A.config.Q(:))),toc(t0));
ix=@(k) find(abs(Gm.x-Gm.trace_x(k))<1e-6);
t1=tic;H=Gm;H.v(:)=Gm.cap_v;H.rho(:)=Gm.cap_rho;H.Q(:)=Gm.cap_Q;[ha,he]=fd_shot(H,ix(50),P,round(P.fd.r0/P.fd.dx));fprintf('shot 0: %.0fs\n',toc(t1));
t1=tic;[a,e]=fd_shot(Gm,ix(j),P);fprintf('shot %d: %.0fs\n',j,toc(t1));
ha=double(single(ha));he=double(single(he));a=double(single(a));e=double(single(e));
tr=process_zo(a(:,1),ha(:,1),he(:,2),Gm.cap_v,Gm.config,P);ref=A.clean(:,j);
fprintf('trace %d: rel diff to archived MATLAB %.2e, corr %.8f\n',j,norm(tr-ref)/norm(ref),corr(tr,ref));
end
