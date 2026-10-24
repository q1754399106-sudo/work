function run_selftest()
% Package checks: unit algebra, and exact reproduction of archived round-3 outputs
% (round3_20261024/results) through the public entry.
here=fileparts(mfilename('fullpath'));root=fileparts(here);addpath(root);addpath(fullfile(root,'core'),fullfile(root,'source'));
maxNumCompThreads(4);data=fullfile(here,'data');
cases={'F4','Q-30','Strong',-30;'V4','Q+0','Strong',0;'F3','Q+0','Moderate',0};
for k=1:size(cases,1)
 f=fullfile(data,[cases{k,1} '.mat']);if ~isfile(f),fprintf('skip %s (data not present)\n',cases{k,1});continue;end
 in=load(f);ref=load(fullfile(data,sprintf('ref_%s_%s_%s.mat',cases{k,1:3})),'gathers','scan');
 ob=in.obs(strcmp({in.obs.name},cases{k,3}));cm=in.config;cm.Q=reshape(cm.Q,1,[])*(1+cases{k,4}/100);
 a=q_compensate(ob.D,cm,'sigma',ob.sigma);b=q_compensate(ob.D,cm,'sigma',ob.sigma,'mode','qscan');
 e1=norm(a.gather-ref.gathers.Standard,'fro')/norm(ref.gathers.Standard,'fro');
 e2=norm(b.gather-ref.gathers.QScan_Standard,'fro')/norm(ref.gathers.QScan_Standard,'fro');
 e3=abs(b.qscan_factor-ref.scan.factor);
 fprintf('%s %s %s: standard diff %.2g, qscan+standard diff %.2g, factor diff %.2g, flags %d/%d\n',cases{k,1:3},e1,e2,e3,a.checks.optimization_flags,b.checks.optimization_flags);
 assert(e1<1e-10&&e2<1e-10&&e3==0,'Package does not reproduce the archived round-3 result.');
end
z=q_compensate(zeros(size(ob.D)),cm,'sigma',ones(1,size(ob.D,2)));assert(all(z.gather(:)==0),'All-dead input must give zero output.');
fprintf('All package self-tests passed.\n');
end
