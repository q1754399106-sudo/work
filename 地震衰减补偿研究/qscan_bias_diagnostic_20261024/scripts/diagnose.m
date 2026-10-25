function [r,x1,xs]=diagnose(P,D,cm,Qused,offset,boot_seed)
% Truth-free diagnostics for "will QScan beat FullCalRefit on this gather?". Inputs are the
% observed gather and the public configuration only.
%  B0: bootstrap world where the used Q is right (truth = FullCalRefit output);
%  G : round-3 gate, relative objective gain over the grid point nearest 1;
%  CV: hold-out prediction of the data (every fifth 20 ms block held out, one fold);
%  RS: residual energy above the noise floor, and residual lag-1 correlation beyond the noise model.
% The mirror world where QScan's Q is right (B1) was dropped for run time: in the development
% probe it was strongly negative in all four scenes, i.e. carried no sign information.
nt=cm.nt;t0=tic;sig=estimate_sigmas(D,cm.t);sel=qscan_select(D,sig,cm,Qused,offset);c=P.grid(sel.best);
x1=full_refit(D,sig,cm,Qused,offset);xs=full_refit(D,sig,cm,Qused*c,offset);
[~,k1]=min(abs(log(P.grid)));gain=(sel.objective(k1)-sel.objective(sel.best))/sel.objective(k1);
[phi,sigb,tend]=noise_model(D,cm);seeds=boot_seed+(1:P.K);tm(1)=toc(t0);t0=tic;
H0=bootstrap_world(x1,Qused,Qused,Qused*c,cm,offset,phi,sigb,seeds);tm(2)=toc(t0);t0=tic;
obs=mod(ceil((1:nt)'/5)-1,5)~=0;cv=[0 0];Qs={Qused,Qused*c};
for q=1:2,Dp=masked_refit(D,sig,cm,Qs{q},offset,obs);cv(q)=sum((D(~obs,:)-Dp(~obs,:)).^2,'all');end
tm(3)=toc(t0);
floor_energy=sum(nt*sigb.^2);sm=cm.t>tend;rs=[0 0];rc=[0 0];X={x1,xs};
for q=1:2
 cc=cm;cc.Q=Qs{q};d=prune_dictionary(cc,offset);R=D-d.G*X{q};rs(q)=sum(R.^2,'all')/floor_energy-1;
 Z=R(sm,:)-mean(R(sm,:),1);rc(q)=sum(Z(2:end,:).*Z(1:end-1,:),'all')/sum(Z.^2,'all')-phi;
end
r=struct('qscan_factor',c,'gate_gain',gain,'phi_hat',phi,'noise_window_end',tend,'sigma_mean',mean(sig), ...
 'B0_rel_delta_mse',H0.rel_delta_mse,'B0_delta_bias2',H0.bias2(2)-H0.bias2(1),'B0_delta_var',H0.var(2)-H0.var(1), ...
 'B0_mse_base',H0.mse(1),'cv_rel_delta',(cv(2)-cv(1))/cv(1),'rs_delta',rs(2)-rs(1),'rc_delta',rc(2)-rc(1), ...
 'seconds_fit',tm(1),'seconds_B0',tm(2),'seconds_cv',tm(3));
end
