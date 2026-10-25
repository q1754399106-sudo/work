function check_calibration()
% Numerical step convergence for the finite-difference calibration Jacobian.
% Does not change or re-fit any sealed prediction.
root=fileparts(fileparts(mfilename('fullpath')));repo=fileparts(root);
addpath(fullfile(repo,'bae_q_20261024','source'));maxNumCompThreads(4);
s=load(fullfile(root,'inputs','complex_Moderate.mat'));p=s.calibration.public;c=s.cm;
p.t=p.t(:);p.vp=reshape(p.vp,1,[]);p.h=reshape(p.h,1,[]);
times=cumsum(2*p.h./p.vp(1:3));mask=false(p.nt,1);for t=times,mask=mask|abs(p.t-t)<=.12;end
rows={};
for factor=unique([1 s.scaninfo.selected_factor])
 q=c.Q*factor;
 for j=1:numel(p.offset)
  ni=estimate_trace_snr(s.calibration.D(:,j),p.t);sig=ni.sigma;
  C=p1_layered_calibration_templates(p,p.offset(j),q,c.fpeak);C=C(mask,:)/sig;
  coef=pinv(C)*(s.calibration.D(mask,j)/sig);[U,S,~]=svd(C,'econ');U=U(:,diag(S)>1e-10*norm(C,2));
  for l=1:5
   ref=column(1e-5);
   for h=[1e-3 1e-4]
    v=column(h);rel=norm(v-ref)/max(norm(ref),1e-12);
    rows(end+1,:)={factor,j,l,h,rel}; %#ok<AGROW>
   end
  end
 end
end
T=cell2table(rows,'VariableNames',{'Q_factor','calibration_trace','parameter','step','relative_to_step_1e5'});
writetable(T,fullfile(root,'analysis','calibration_fd_convergence.csv'));
assert(max(T.relative_to_step_1e5(T.step==1e-4))<1e-5);
fprintf('Calibration FD stability max %.6g\n',max(T.relative_to_step_1e5(T.step==1e-4)));
 function v=column(h)
  qp=q;qm=q;fp=c.fpeak;fm=c.fpeak;
  if l<=4,qp(l)=qp(l)*exp(h);qm(l)=qm(l)*exp(-h);else,fp=fp*exp(h);fm=fm*exp(-h);end
  P=p1_layered_calibration_templates(p,p.offset(j),qp,fp);M=p1_layered_calibration_templates(p,p.offset(j),qm,fm);
  v=(P(mask,:)-M(mask,:))/(2*h)*coef/sig;v=v-U*(U'*v);
 end
end
