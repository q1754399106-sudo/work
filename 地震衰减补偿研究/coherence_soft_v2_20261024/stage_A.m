function stage_A()
% Registered stage A: extended lambda grid and hybrid deletion of isolated groups.
ref=cand('v1_1_8',1,8,0);c=ref([]);
for lmin=[.5 1],for lmax=[8 12 16 32],if lmin==1&&lmax==8,continue;end,c(end+1)=cand(sprintf('lam_%g_%g',lmin,lmax),lmin,lmax,0);end,end %#ok<AGROW>
for tau=[.1 .2],for lmax=[8 16],c(end+1)=cand(sprintf('del%g_1_%g',tau,lmax),1,lmax,tau);end,end %#ok<AGROW>
run_stage('A',c,ref);
end

function x=cand(label,lmin,lmax,tau)
x=struct('label',label,'opts',struct('lambda_min',lmin,'lambda_max',lmax,'tau',tau,'kappa',0),'cfun',[]);
end
