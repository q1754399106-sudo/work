function bench_fd()
% Timing of fd2d_visco on candidate grids (no physics checked here).
maxNumCompThreads(4);
for g=[4 5 4 5]
 nz=round(2950/g)+40;nx=round(1000/g)+40;M.v=2500*ones(nz,nx);M.v(round(nz/2):end,:)=3300;M.rho=.31*M.v.^.25;M.Q=60*ones(nz,nx);
 dt=g*1e-4;opt=struct('dx',g,'dt',dt,'nt',200,'fref',100,'npml',20,'src',[21 round(nx/2)],'rec',[21 round(nx/2)]);
 q=zeros(200,1);q(1:20)=1;
 o2=opt;o2.nt=20;fd2d_visco(M,q,o2);   % warm-up
 tic;fd2d_visco(M,q,opt);tv=toc;M.Q(:)=Inf;tic;fd2d_visco(M,q,opt);te=toc;
 steps=ceil(1.85/dt);
 fprintf('dx=%g m: grid %dx%d (%d cells), dt=%.2g ms, %d steps for 1.85 s | visco %.1f ms/step (%.0f s/run), elastic %.1f ms/step (%.0f s/run)\n', ...
  g,nz,nx,nz*nx,dt*1e3,steps,tv/200*1e3,tv/200*steps,te/200*1e3,te/200*steps);
end
end
