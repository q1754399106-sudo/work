function [r,x1,xs]=diagnose2(P,D,cm,Qused,offset,boot_seed)
% Truth-free B0 with three choices of the bootstrap world's truth (same noise draws in all):
%  B0  : FullCalRefit output x1 (qscan_bias_diagnostic_20261024, reproduced exactly);
%  B0M : bootstrap mean m* = average FullCalRefit output over the B0 world's draws;
%  B0C : coherence-shrunk output (round-3 Standard) on the observed gather.
% Also the round-3 gate G and the variance-only call V of each world.
t0=tic;sig=estimate_sigmas(D,cm.t);sel=qscan_select(D,sig,cm,Qused,offset);c=P.grid(sel.best);
x1=full_refit(D,sig,cm,Qused,offset);xs=full_refit(D,sig,cm,Qused*c,offset);
[~,k1]=min(abs(log(P.grid)));gain=(sel.objective(k1)-sel.objective(sel.best))/sel.objective(k1);
[phi,sigb,tend]=noise_model(D,cm);seeds=boot_seed+(1:P.K);
[XC,xcheck,ci]=coherence_reference(P,D,sig,cm,Qused,offset);assert(isequal(xcheck,x1),'coherence path changed FullCalRefit.');
H0=bootstrap_world2(x1,Qused,Qused,Qused*c,cm,offset,phi,sigb,seeds);
HM=bootstrap_world2(H0.mean_a,Qused,Qused,Qused*c,cm,offset,phi,sigb,seeds);
HC=bootstrap_world2(XC,Qused,Qused,Qused*c,cm,offset,phi,sigb,seeds);
e=@(X)sum(X(:).^2);
r=struct('qscan_factor',c,'gate_gain',gain,'phi_hat',phi,'noise_window_end',tend, ...
 'B0_rel_delta_mse',H0.rel_delta_mse,'B0_delta_bias2',H0.bias2(2)-H0.bias2(1),'B0_delta_var',H0.var(2)-H0.var(1), ...
 'B0M_rel_delta_mse',HM.rel_delta_mse,'B0M_delta_bias2',HM.bias2(2)-HM.bias2(1),'B0M_delta_var',HM.var(2)-HM.var(1), ...
 'B0C_rel_delta_mse',HC.rel_delta_mse,'B0C_delta_bias2',HC.bias2(2)-HC.bias2(1),'B0C_delta_var',HC.var(2)-HC.var(1), ...
 'energy_mean_over_x1',e(H0.mean_a)/e(x1),'energy_coh_over_x1',e(XC)/e(x1),'mean_coherence',ci.mean_coherence, ...
 'seconds',toc(t0));
end
