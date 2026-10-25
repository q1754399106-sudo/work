function sc=load_obs(P,set,model,condition)
% Solver-side inputs only. set: 'bae' = bae_q F1/F2, 'qc' = qscan_criterion F3/F4,
% 'r3' = round3 F3/F4 (two different geologies share the names F3/F4).
switch set
 case 'bae',in=fullfile(P.in,'bae_q_20261024','inputs','independent',[model '.mat']);
 case 'qc',in=fullfile(P.in,'qscan_criterion_20261024','inputs','confirm',[model '.mat']);
 case 'r3',in=fullfile(P.in,'round3_20261024','dev',[model '.mat']);
end
s=load(in);assert(~isfield(s,'tr'));cm=s.config;cm.t=cm.t(:);k=strcmp({s.obs.name},condition);assert(nnz(k)==1);
sc=struct('set',set,'model',model,'condition',condition,'cm',cm,'D',s.obs(k).D, ...
 'sigma',reshape(s.obs(k).sigma,1,[]),'Q0',reshape(cm.Q,1,[]),'offset',cm.solver_offset);
end
