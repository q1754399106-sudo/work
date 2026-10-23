function p1_regression_checks(outdir)
% Targeted regressions for the observed bugs, independent of acceptance gains.
root=fileparts(fileparts(mfilename('fullpath')));
spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json'))); p1=spec.p1_hyperparameters;
cm=config_model(); cm.nt=161; cm.t=(0:160)'*cm.dt; cm.nfft=1024;
cm.Q=40*ones(1,4); cm.fpeak=55;
cm.log_parameter_covariance=[.04,-.006;-.006,.01];
[nodes,w]=p1_parameter_nodes(cm); logs=log(nodes); meanlog=logs*w';
dc=logs-meanlog; recovered=(dc.*w)*dc';
assert(norm(meanlog-log([40;55]))<1e-12 && norm(recovered-cm.log_parameter_covariance,'fro')<1e-12);
n=cm.nt; model=struct('G',eye(n),'F',eye(n),'K',eye(n),'gram_direct',eye(n),'gram_split',eye(n), ...
 'direct',struct('U',eye(n),'V',eye(n),'z',ones(n,1)), ...
 'split',struct('U',eye(n),'V',eye(n),'z',ones(n,1)));
models=repmat({model},1,11); R=sparse(n*11,n*11);
rng(58121,'twister'); D=1e-3*randn(n,11); D(:,end)=D(:,end)*.01;
[a,prior]=p1_solve_branch(D,models,ones(1,11),cm,'Section',R,'CombinedUncertainty_NoStructure');
assert(nnz(prior.extra_diag)>0,'Regression fixture must activate tail regularization');
assert(a.direct_record.normal_residual<=2e-7 && a.split_record.normal_residual<=2e-7);
[~,qprior]=p1_solve_branch(D,models,ones(1,11),cm,'Section',R,'QModelErrorOnly');
assert(nnz(qprior.extra_diag)==0 && nnz(qprior.tail_mask)==0,'Pure Q branch leaked noise protection');
[r_unc,~,~,~,~]=p1_q_model_error(D,models,cm,'Section',true,false,p1);
cz=cm; cz.log_parameter_covariance=zeros(2);
[r_zero,~,~,~,~]=p1_q_model_error(D,models,cz,'Section',true,false,p1);
f=(0:floor(n/2))'/(n*cm.dt); Vf=abs(fft(model.direct.V,n,1)); Vf=Vf(1:numel(f),:).^2;
fm=((f'*Vf)./max(sum(Vf,1),realmin))';
legacy=p1_modal_penalty(ones(n,1),fm,40,55,p1);
assert(norm(r_zero{1}-legacy)<1e-12,'Zero uncertainty must recover point penalty');
weight_change=norm(r_unc{1}-r_zero{1}); assert(weight_change>1e-6,'Uncertainty is not used');
small=config_model(); small.nt=161; small.t=(0:160)'*small.dt; small.nfft=1024;
pub=uq_public_config(small,55);
est=struct('effective_Q',40,'estimated_frequency',55,'log_parameter_covariance',zeros(2));
ca=uq_inverse_config(pub,est); ma=prepare_models_v2(750,ca);
est.effective_Q=80; cb=uq_inverse_config(pub,est); mb=prepare_models_v2(750,cb);
operator_change=norm(ma{1}.G-mb{1}.G,'fro'); assert(operator_change>1e-4);
checks=struct('pass',true,'cubature_covariance_error',norm(recovered-cm.log_parameter_covariance,'fro'), ...
 'active_tail_diagonal_entries',nnz(prior.extra_diag),'no_structure_normal_residual',a.direct_record.normal_residual, ...
 'pure_q_noise_entries',nnz(qprior.extra_diag),'uncertainty_weight_change',weight_change, ...
 'zero_uncertainty_point_penalty_error',norm(r_zero{1}-legacy),'recalibrated_operator_change',operator_change);
uq_writejson(fullfile(outdir,'regression_checks.json'),checks);
fprintf('PASS: repaired no-structure solve, pure-Q isolation, uncertainty propagation, operator refresh\n');
end
