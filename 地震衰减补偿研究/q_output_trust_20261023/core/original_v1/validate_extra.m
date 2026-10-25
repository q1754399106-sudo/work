function validate_extra()
c=config_model(); load(fullfile(c.output,'experiment.mat'),'D','names','c');
E=event_metrics(D,names,c); writetable(E,fullfile(c.output,'event_metrics.csv'));
% Independent constant-Q, homogeneous-medium analytic check.
cc=c; cc.vp(:)=2200; cc.Q(:)=60;
t=[0.3;0.7;1.0]; x=900;
[a,b]=ray_absorption(t,x,cc);
analytic=sqrt(t.^2+(x/2200)^2)/60;
assert(norm(a+b-analytic)<1e-12,'Homogeneous ray mismatch');
% For constant absorption, a delayed delta should match its padded impulse.
cc.nt=151; cc.nfft=2048;
G=attenuation_matrix(0.01*ones(cc.nt,1),cc);
assert(norm(G(2:end,2:end)-G(1:end-1,1:end-1),'fro')<1e-10);
fid=fopen(fullfile(c.output,'checks.txt'),'a');
fprintf(fid,'PASS: homogeneous analytic ray loss and constant-loss shift invariance.\n'); fclose(fid);
fprintf('Extra checks passed.\n');
end
