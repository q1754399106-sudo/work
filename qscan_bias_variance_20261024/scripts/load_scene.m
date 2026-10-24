function sc=load_scene(P,s)
% Solver inputs exactly as list_scenes('ind'); D0 (noise-free gather) and the target are
% loaded for data synthesis and scoring only and are never passed to a solver.
spec=P.scenes(s);in=load(fullfile(P.bae,'inputs','independent',[spec.model '.mat']));
cm=in.config;cm.t=cm.t(:);k=find(strcmp({in.obs.name},spec.condition));
sc=spec;sc.cm=cm;sc.D=in.obs(k).D;sc.sigma=reshape(in.obs(k).sigma,1,[]);
sc.Q0=reshape(cm.Q,1,[]);sc.offset=cm.solver_offset;
t=load(fullfile(P.bae,'truth','independent',[spec.model '_truth.mat']));
sc.T=t.tr(strcmp({t.tr.name},spec.condition)).target;sc.D0=t.clean;
end
