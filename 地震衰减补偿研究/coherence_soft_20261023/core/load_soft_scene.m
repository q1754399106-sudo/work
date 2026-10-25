function scene=load_soft_scene(softroot,which,errorpct,name,D)
% Assemble the whole-gather FullCalRefit state saved by run_soft_trial.
nx=size(D,2);R0=zeros(size(D));supports=cell(1,nx);E=supports;V=supports;S0=R0;soft=R0;
for block=1:nx/11
 ids=(block-1)*11+(1:11);
 c=load(fullfile(softroot,'cache',sprintf('%s_Q%+d_%s_%02d.mat',which,errorpct,name,block)),'initial','out');
 R0(:,ids)=c.initial.reflectivity;S0(:,ids)=c.initial.gather;soft(:,ids)=c.out.gather;
 supports(ids)=c.initial.support;E(ids)=c.out.group_energy;V(ids)=c.out.group_noise_energy;
end
scene=struct('R0',R0,'S0',S0,'soft',soft,'supports',{supports},'E',{E},'V',{V});
end
