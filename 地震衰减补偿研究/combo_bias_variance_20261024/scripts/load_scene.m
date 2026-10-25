function sc=load_scene(P,s)
% Public config of the scene's model, plus the noise-free gather D0 and the target. D0 is used
% only to synthesise observations and the target only to score; neither reaches a solver.
sc=P.scenes(s);in=load(fullfile(P.bae,'inputs','independent',[sc.model '.mat']),'config');
cm=in.config;cm.t=cm.t(:);sc.cm=cm;sc.Q0=reshape(cm.Q,1,[]);sc.Qused=sc.Q0*sc.q_scale;sc.offset=cm.solver_offset;
t=load(fullfile(P.bae,'truth','independent',[sc.model '_truth.mat']),'tr','clean');
assert(isequal(t.tr(1).target,t.tr(2).target),'Targets differ between conditions.');
sc.T=t.tr(1).target;sc.D0=t.clean;
end
