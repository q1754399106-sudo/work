function score_dev()
% Development scoring and the pre-registered selection rule (PLAN_ZH.md).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));
opts=jsondecode(fileread(fullfile(root,'SPEC.json')));errs=reshape(opts.errors_percent,1,[]);soft=fullfile(fileparts(root),'soft_refit_20261023');
rows={};flags=0;normal=0;
for which={'primary','complex'}
 [cm,ref,input]=load_dev(soft,which{1});
 for k=1:numel(ref.datasets)
  data=ref.datasets(k);live=any(data.D~=0,1);T=data.target;
  for e=errs
   L=load(fullfile(root,'results','dev',sprintf('%s_Q%+d_%s.mat',which{1},e,data.name)));
   flags=flags+L.diag.optimization_flags;normal=max(normal,L.diag.max_normal_residual);
   for m=reshape(fieldnames(L.g),1,[])
    A=double(L.g.(m{1}));rel=norm(A(:,live)-T(:,live),'fro')/norm(T(:,live),'fro');weak=NaN;
    if strcmp(which{1},'primary'),weak=score_weak(A,input.reflectivity,input.event_indices_zero_based+1,[2 4 6 8],cm.fpeak,cm.dt);end
    rows(end+1,:)={which{1},data.name,e,m{1},rel,weak}; %#ok<AGROW>
   end
  end
 end
end
M=cell2table(rows,'VariableNames',{'model','condition','Q_percent','method','relative_error','weak_error'});
writetable(M,fullfile(root,'analysis','dev_metrics.csv'));
methods=unique(M.method);combos=unique(M(:,{'model','condition'}));S=table();
for q=1:numel(methods)
 for c=1:height(combos)
  sel=strcmp(M.model,combos.model{c})&strcmp(M.condition,combos.condition{c});
  base=M(sel&strcmp(M.method,'FullCalRefit'),:);mine=M(sel&strcmp(M.method,methods{q}),:);
  wc=max(mine.relative_error)/max(base.relative_error);q0=mine.relative_error(mine.Q_percent==0)/base.relative_error(base.Q_percent==0);
  w0=mine.weak_error(mine.Q_percent==0)/base.weak_error(base.Q_percent==0);
  S=[S;table(methods(q),combos.model(c),combos.condition(c),wc,q0,w0,'VariableNames',{'method','model','condition','wc_ratio','q0_ratio','weak_q0_ratio'})]; %#ok<AGROW>
 end
end
writetable(S,fullfile(root,'analysis','dev_ratios.csv'));
selected=struct();
for family={'BAE','BAECovOnly','QEnsembleMean'}
 best='';bestwc=Inf;cands=struct('name',{},'mean_wc',{},'mean_q0',{},'mean_weak_q0',{},'feasible',{});
 for s=reshape(opts.s_grid,1,[])
  name=sprintf('%s_s%03d',family{1},round(100*s));r=S(strcmp(S.method,name),:);
  mwc=mean(r.wc_ratio);mq0=mean(r.q0_ratio);mw=mean(r.weak_q0_ratio(strcmp(r.model,'primary')));
  ok=mq0<=1.02&&mw<=1.10;cands(end+1)=struct('name',name,'mean_wc',mwc,'mean_q0',mq0,'mean_weak_q0',mw,'feasible',ok); %#ok<AGROW>
  if ok&&mwc<bestwc,best=name;bestwc=mwc;end
 end
 selected.(family{1})=struct('selected',best,'candidates',cands);
end
selected.optimization_flags=flags;selected.max_normal_residual=normal;
uq_writejson(fullfile(root,'analysis','dev_selection.json'),selected);
disp(S);disp(jsonencode(selected));
end
