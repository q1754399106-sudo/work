function score_crit(split)
% Truth-side scoring of every registered criterion (PLAN_ZH.md).
root=fileparts(mfilename('fullpath'));bae=fullfile(fileparts(root),'bae_q_20261024');
addpath(fullfile(root,'core'),fullfile(bae,'core'),fullfile(bae,'source'));spec=jsondecode(fileread(fullfile(root,'SPEC.json')));
crit=spec.criteria;files=dir(fullfile(root,'results',split,'*.mat'));rows={};
for f=reshape(files,1,[])
 tok=regexp(f.name,'^(.*)_Q([+-]\d+)_(.*)\.mat$','tokens','once');model=tok{1};e=str2double(tok{2});name=tok{3};
 L=load(fullfile(f.folder,f.name));st=L.st;g=st.grid;base=find(abs(g-1)<1e-12);
 if strcmp(split,'field')
  for c=1:numel(crit),k=select_factor(st,crit{c});rows(end+1,:)={model,name,e,crit{c},g(k),st.used_Q(2)*g(k),NaN,NaN,NaN,NaN,k==1||k==numel(g)};end %#ok<AGROW>
  continue
 end
 T=truth(root,bae,split,model,name);live=st.live;rel=zeros(numel(g),1);
 for q=1:numel(g),A=double(st.gathers(:,:,q));rel(q)=norm(A(:,live)-T(:,live),'fro')/norm(T(:,live),'fro');end
 [ro,ko]=min(rel);
 for c=1:numel(crit)
  k=select_factor(st,crit{c});
  rows(end+1,:)={model,name,e,crit{c},g(k),NaN,abs(log(g(k)/g(ko))),rel(k),ro,rel(base),k==1||k==numel(g)}; %#ok<AGROW>
 end
end
M=cell2table(rows,'VariableNames',{'model','condition','Q_percent','criterion','factor','implied_Q2','log_factor_error','rel_error','oracle_rel_error','base_rel_error','at_grid_edge'});
writetable(M,fullfile(root,'analysis',[split '_criteria.csv']));if strcmp(split,'field'),disp(M);return;end
S=struct();
for c=1:numel(crit)
 r=M(strcmp(M.criterion,crit{c}),:);combos=unique(r(:,{'model','condition'}));wc=zeros(height(combos),1);q0=wc;
 for i=1:height(combos)
  s=r(strcmp(r.model,combos.model{i})&strcmp(r.condition,combos.condition{i}),:);
  wc(i)=max(s.rel_error)/max(s.base_rel_error);q0(i)=s.rel_error(s.Q_percent==0)/s.base_rel_error(s.Q_percent==0);
 end
 lay=ismember(combos.model,{'primary','complex','V1','V2','V3'});
 S.(crit{c})=struct('mean_log_factor_error',mean(r.log_factor_error),'mean_ratio_to_oracle',mean(r.rel_error./r.oracle_rel_error), ...
  'mean_wc_ratio',mean(wc),'mean_q0_ratio',mean(q0),'mean_wc_ratio_layered',mean(wc(lay)),'mean_q0_ratio_layered',mean(q0(lay)), ...
  'edge_hits',sum(r.at_grid_edge),'cases',height(r));
end
uq_writejson(fullfile(root,'analysis',[split '_summary.json']),S);disp(jsonencode(S));
end

function T=truth(root,bae,split,model,name)
switch model
 case {'primary','complex'}
  [~,ref]=load_dev(fullfile(fileparts(root),'soft_refit_20261023'),model);T=ref.datasets(strcmp({ref.datasets.name},name)).target;
 case {'F1','F2'},t=load(fullfile(bae,'truth','independent',[model '_truth.mat']));T=t.tr(strcmp({t.tr.name},name)).target;
 otherwise,t=load(fullfile(root,'truth',split,[model '_truth.mat']));T=t.tr(strcmp({t.tr.name},name)).target;
end
end
