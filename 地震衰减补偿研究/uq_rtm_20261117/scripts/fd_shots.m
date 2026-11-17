function fd_shots(name,ks)
% Viscoacoustic (GSLS) and elastic (Q = Inf) shots of model name for shot indices ks; resumable.
% Shot k: source at (P.izs, P.ixs(k)), receivers at (P.izs, P.ixr). Output cache/shots_<name>/sNNN.mat:
% att, ela (P.nt x 300 single, t = (1:P.nt)*P.dt, low-passed below 100 Hz before decimation).
% name 'T0' (homogeneous top layer): one shot at x = 8 m recorded at every interior column (direct wave).
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model(name,P);
out=fullfile(P.cache,['shots_' name]);if ~isfolder(out),mkdir(out);end
q=ricker(P);opt=struct('dx',P.dx,'dt',P.fd.dt,'nt',P.fd.nt,'fref',P.fd.fref,'npml',P.npml);
if strcmp(name,'T0'),ks=0;end
for k=ks
 fn=fullfile(out,sprintf('s%03d.mat',k));if isfile(fn),continue;end
 if k==0
  opt.src=[P.izs P.npml+3];cols=P.in{2};
 else
  opt.src=[P.izs P.ixs(k)];cols=P.ixr;
 end
 opt.rec=[P.izs*ones(numel(cols),1) cols(:)];
 t0=tic;att=decim(fd2d_visco(M,q,opt));Me=M;Me.Q(:)=Inf;ela=decim(fd2d_visco(Me,q,opt));
 save(fn,'att','ela');fprintf('%s shot %d: %.0f s\n',name,k,toc(t0));
end
end

function q=ricker(P)
t=(1:P.fd.nt)'*P.fd.dt-P.fd.t0;a=(pi*P.fd.fpeak*t).^2;q=(1-2*a).*exp(-a);
end

function y=decim(x)
n=size(x,1);f=(0:n-1)'/(n*4e-4);f=min(f,1/4e-4-f);
w=ones(n,1);b=f>80&f<110;w(b)=.5*(1+cos(pi*(f(b)-80)/30));w(f>=110)=0;
y=real(ifft(fft(x).*w));y=single(y(10:10:end,:));
end
