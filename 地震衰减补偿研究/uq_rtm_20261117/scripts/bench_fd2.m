function bench_fd2()
maxNumCompThreads(4);
nz=440;nx=640;M.v=2500*ones(nz,nx);M.v(200:end,:)=3300;M.rho=2000*ones(nz,nx);M.Q=60*ones(nz,nx);
rec=[22*ones(300,1) round(linspace(30,610,300))'];
opt=struct('dx',4,'dt',4e-4,'nt',200,'fref',100,'npml',25,'src',[22 320],'rec',rec);q=zeros(200,1);q(1:20)=1;
o2=opt;o2.nt=20;fd2d_visco(M,q,o2);
tic;fd2d_visco(M,q,opt);tv=toc;M.Q(:)=Inf;tic;fd2d_visco(M,q,opt);te=toc;
fprintf('visco %.1f ms/step -> %.0f s for 5000 steps; elastic %.1f ms/step -> %.0f s\n',tv/200*1e3,tv/200*5000,te/200*1e3,te/200*5000);
end
