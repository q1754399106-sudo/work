function diag_datalayer(sets)
% Data-only layered QScan check (rule in analysis/ORACLE_AND_DECISION_RULE.md).
% Factors chosen by the lasso objective J only; truth used afterwards for error.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(2);
grid=exp(linspace(log(.6),log(2.2),13));traces=1:3:99;rows={};started=tic;
for v=sets
 in=load(fullfile(root,'dev',[v{1} '.mat']));T=load(fullfile(root,'dev',[v{1} '_truth.mat']),'tr');
 cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);offset=450;if isfield(cm,'solver_offset'),offset=cm.solver_offset;end
 D=in.obs(1).D(:,traces);sig=reshape(in.obs(1).sigma(traces),1,[]);tg=T.tr(1).target(:,traces);
 for e=[-30 0 30]
  used=cm;used.Q=Q0*(1+e/100);obj=@(c)lasso_objective(D,sig,used,c,offset);
  g=arrayfun(@(x)obj(x*ones(1,4)),grid);[~,ib]=min(g);c=grid(ib)*ones(1,4);best=g(ib);
  for sweep=1:2
   changed=false;
   for L=1:4
    for x=grid
     if x==c(L),continue;end
     trial=c;trial(L)=x;j=obj(trial);if j<best,best=j;c=trial;changed=true;end
    end
   end
   if ~changed,break;end
  end
  eg=fullcal_error(D,sig,used,grid(ib)*ones(1,4),offset,tg);el=fullcal_error(D,sig,used,c,offset,tg);
  rows{end+1}=struct('set',v{1},'Q_percent',e,'data_global_factor',grid(ib),'data_layer_factors',c,'data_global_error',eg,'data_layered_error',el,'ratio',el/eg); %#ok<AGROW>
  fprintf('%s Q%+d: data global %.4f (c=%.2f) | data layered %.4f (c=%s) | ratio %.3f | %.0fs\n',v{1},e,eg,grid(ib),el,mat2str(round(c,2)),el/eg,toc(started));
 end
end
uq_writejson(fullfile(root,'analysis',sprintf('datalayer_%s.json',strjoin(sets,'_'))),[rows{:}]);
end

function J=lasso_objective(D,sigma,used,c,offset)
% Same objective as qscan_refit.m (gather-summed, sigma-normalised).
n=size(D,1);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;cm=used;cm.Q=used.Q.*c;d=prune_dictionary(cm,offset);
X=refit_lasso(d.Gram,d.B'*D(:,live),lambda(live),d.L);res=D(:,live)-d.B*X;
J=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);
end

function e=fullcal_error(D,sig,used,c,offset,tg)
cm=used;cm.Q=used.Q.*c;d=prune_dictionary(cm,offset);nx=size(D,2);S=zeros(size(D));
for b=1:11:nx,ids=b:min(b+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sig(ids));S(:,ids)=o.gather;end
e=norm(S-tg,'fro')/norm(tg,'fro');
end
