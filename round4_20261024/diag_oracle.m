function diag_oracle(sets)
% Upper-bound diagnosis (truth used to choose factors; diagnosis only).
% FullCalRefit error with: nominal Q, best common factor (oracle global),
% best per-layer factors (oracle layered, coordinate descent from the global
% optimum). Every third trace (33 traces) to bound cost; per-trace solves are
% independent, so the subset changes no individual solution.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(2);
grid=exp(linspace(log(.6),log(2.2),13));traces=1:3:99;rows={};started=tic;
for v=sets
 in=load(fullfile(root,'dev',[v{1} '.mat']));T=load(fullfile(root,'dev',[v{1} '_truth.mat']),'tr');
 cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);offset=450;if isfield(cm,'solver_offset'),offset=cm.solver_offset;end
 k=1;D=in.obs(k).D(:,traces);sig=reshape(in.obs(k).sigma(traces),1,[]);tg=T.tr(k).target(:,traces);
 for e=[-30 0 30]
  used=cm;used.Q=Q0*(1+e/100);err=@(c)fullcal_error(D,sig,used,c,offset,tg);
  nominal=err(ones(1,4));g=arrayfun(@(x)err(x*ones(1,4)),grid);[gbest,ib]=min(g);
  c=grid(ib)*ones(1,4);best=gbest;
  for sweep=1:2
   before=best;
   for L=1:4
    for x=grid
     if x==c(L),continue;end
     trial=c;trial(L)=x;v_=err(trial);if v_<best,best=v_;c=trial;end
    end
   end
   if best>.99*before,break;end
  end
  rows{end+1}=struct('set',v{1},'Q_percent',e,'nominal',nominal,'oracle_global',gbest,'global_factor',grid(ib), ...
   'oracle_layered',best,'layer_factors',c,'layered_vs_global',best/gbest,'global_vs_nominal',gbest/nominal); %#ok<AGROW>
  fprintf('%s Q%+d: nominal %.4f | global %.4f (c=%.2f) | layered %.4f (c=%s) | layered/global %.3f | %.0fs\n', ...
   v{1},e,nominal,gbest,grid(ib),best,mat2str(round(c,2)),best/gbest,toc(started));
 end
end
uq_writejson(fullfile(root,'analysis',sprintf('oracle_%s.json',strjoin(sets,'_'))),[rows{:}]);
end

function e=fullcal_error(D,sig,used,c,offset,tg)
cm=used;cm.Q=used.Q.*c;d=prune_dictionary(cm,offset);nx=size(D,2);S=zeros(size(D));
for b=1:11:nx
 ids=b:min(b+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sig(ids));S(:,ids)=o.gather;
end
e=norm(S-tg,'fro')/norm(tg,'fro');
end
