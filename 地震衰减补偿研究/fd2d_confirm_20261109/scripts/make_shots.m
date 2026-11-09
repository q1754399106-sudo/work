function make_shots(m,first,last)
% Worker: zero-offset shots first..last of geology m (traces 1..99). Shot 0 is the laterally
% homogeneous cap model (direct-wave reference at the source node, signature reference r0 below
% it). Raw source-node pressure of the attenuated and the elastic run is kept per shot in
% cache/shots/<model>/; existing shots are skipped, so the worker can be restarted.
P=fd_setup();maxNumCompThreads(P.threads);sp=P.models(m);Gm=fd_geology(sp,P);
dir_=fullfile(P.root,'cache','shots',sp.name);if ~isfolder(dir_),mkdir(dir_);end
ix=@(j) find(abs(Gm.x-Gm.trace_x(j))<1e-6);started=tic;
for j=first:last
 out=fullfile(dir_,sprintf('s%03d.mat',j));if isfile(out),continue;end
 t0=tic;
 if j==0
  H=Gm;H.v(:)=Gm.cap_v;H.rho(:)=Gm.cap_rho;H.Q(:)=Gm.cap_Q;[a,e]=fd_shot(H,ix(50),P,round(P.fd.r0/P.fd.dx));
 else
  [a,e]=fd_shot(Gm,ix(j),P);
 end
 a=single(a);e=single(e);save([out '.part'],'a','e','-v7');movefile([out '.part'],out);
 fprintf('%s shot %d (%.0fs, total %.0fs)\n',sp.name,j,toc(t0),toc(started));
end
fprintf('DONE %s %d-%d %.0fs\n',sp.name,first,last,toc(started));
end
