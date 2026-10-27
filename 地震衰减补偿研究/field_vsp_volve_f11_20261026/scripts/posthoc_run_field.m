function posthoc_run_field(varargin)
% POST-HOC (not registered; RESULTS_ZH.md): same code with Q from analysis/posthoc_vsp_q.json and separate outputs.
% Runs of PLAN_ZH.md ("方法"): combo_outputs on the prepared curtain with Q_VSP x (1+e/100).
% Arguments are run ids (V0, Vm50, Vm30, Vp30, Vp50); a run whose result exists is skipped.
P=setup_paths();maxNumCompThreads(P.threads);root=P.root;
in=load(fullfile(root,'inputs','posthoc_curtain.mat'));vq=jsondecode(fileread(fullfile(root,'analysis','posthoc_vsp_q.json')));
assert(vq.V_executable,'Q_VSP not executable (both segments invalid)');
out=fullfile(root,'results_posthoc');if ~isfolder(out),mkdir(out);end;started=tic;
for id=varargin
 r=P.runs(strcmp({P.runs.id},id{1}));assert(numel(r)==1,'unknown run %s',id{1});
 file=fullfile(out,[r.id '.mat']);if isfile(file),continue;end
 Q=reshape(vq.Q_models.VSP,1,[])*(1+r.e/100);
 [G,info]=combo_outputs(P,in.D,in.sigma,in.cm,Q,0);G.Input=in.D;
 save([file '.part'],'G','info','Q','r','-v7');movefile([file '.part'],file);
 fprintf('%s Q=[%s] QScan factor %.4f flags %d (%.0fs, total %.0fs)\n',r.id,num2str(Q,'%.1f '),info.qscan_factor, ...
  info.optimization_flags,info.seconds,toc(started));
end
fprintf('DONE %s\n',strjoin(varargin,','));
end
