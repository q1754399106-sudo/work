function score_layer(split)
% Truth-side scoring, selection rule (dev) and gates (confirm); field is descriptive.
root=fileparts(mfilename('fullpath'));up=fullfile(root,'upstream');addpath(fullfile(root,'core'),fullfile(up,'bae_core'),fullfile(up,'bae_source'));
files=dir(fullfile(root,'results',split,'*.mat'));rows={};frows={};
if strcmp(split,'dev'),ref=readtable(fullfile(up,'dev_criteria.csv'),'TextType','char');ref=ref(strcmp(ref.criterion,'A0'),:);end
repro=[];
for f=reshape(files,1,[])
 tok=regexp(f.name,'^(.*)_Q([+-]\d+)_(.*)\.mat$','tokens','once');model=tok{1};e=str2double(tok{2});name=tok{3};r=load(fullfile(f.folder,f.name));r=r.r;
 for v=reshape(r.variants,1,[]),frows(end+1,:)={model,name,e,v.name,num2str(v.factors,'%.3f '),num2str(r.used_Q.*v.factors,'%.1f ')};end %#ok<AGROW>
 if strcmp(split,'field'),continue;end
 if strcmp(split,'dev')
  a=ref(strcmp(ref.model,model)&strcmp(ref.condition,name)&ref.Q_percent==e,:);repro(end+1)=abs(a.factor-r.c0)<1e-9; %#ok<AGROW>
 end
 T=truth(root,up,model,name);live=any(T~=0,1)&any(r.gathers.FullCalRefit~=0,1);
 for m=reshape(fieldnames(r.gathers),1,[])
  A=double(r.gathers.(m{1}));rows(end+1,:)={model,name,e,m{1},norm(A(:,live)-T(:,live),'fro')/norm(T(:,live),'fro')}; %#ok<AGROW>
 end
end
F=cell2table(frows,'VariableNames',{'model','condition','Q_percent','variant','factors','applied_Q'});writetable(F,fullfile(root,'analysis',[split '_factors.csv']));
if strcmp(split,'field'),disp(F);return;end
M=cell2table(rows,'VariableNames',{'model','condition','Q_percent','method','rel_error'});writetable(M,fullfile(root,'analysis',[split '_metrics.csv']));
methods=setdiff(unique(M.method),{'FullCalRefit'});S=struct();
for q=1:numel(methods)
 m=methods{q};X=M(strcmp(M.method,m),:);L0=M(strcmp(M.method,'L0'),:);B=M(strcmp(M.method,'FullCalRefit'),:);
 [~,ia,ib]=intersect(strcat(X.model,X.condition,string(X.Q_percent)),strcat(L0.model,L0.condition,string(L0.Q_percent)));ratio=X.rel_error(ia)./L0.rel_error(ib);
 combos=unique(X(:,{'model','condition'}));wc=zeros(height(combos),1);q0=wc;
 for i=1:height(combos)
  s=@(T)T(strcmp(T.model,combos.model{i})&strcmp(T.condition,combos.condition{i}),:);x=s(X);b=s(B);
  wc(i)=max(x.rel_error)/max(b.rel_error);q0(i)=x.rel_error(x.Q_percent==0)/b.rel_error(b.Q_percent==0);
 end
 lay=ismember(combos.model,{'primary','complex'});
 S.(m)=struct('mean_ratio_vs_L0',mean(ratio),'better_than_L0',mean(ratio<1),'mean_wc_ratio',mean(wc),'mean_q0_ratio',mean(q0), ...
  'wc_layered',mean(wc(lay)),'q0_layered',mean(q0(lay)),'cases',numel(ratio));
end
if strcmp(split,'dev')
 S.reproduction_L0_matches_criterion_A0=sprintf('%d/%d',sum(repro),numel(repro));best='';bv=Inf;
 for q=1:numel(methods)
  m=methods{q};if strcmp(m,'L0'),continue;end;s=S.(m);ok=s.wc_layered<=.95&&s.q0_layered<=1.02;
  if ok&&s.mean_ratio_vs_L0<bv,best=m;bv=s.mean_ratio_vs_L0;end
 end
 S.selected=best;S.promoted=~isempty(best)&&bv<.97;
else
 sel=jsondecode(fileread(fullfile(root,'analysis','dev_summary.json')));
 if sel.promoted
  s=S.(sel.selected);S.gates=struct('variant',sel.selected,'H1',s.mean_ratio_vs_L0<=.97,'H2',s.better_than_L0>=.6,'H3',s.mean_wc_ratio<=.95,'H4',s.mean_q0_ratio<=1.02);
  S.gates.all_pass=S.gates.H1&&S.gates.H2&&S.gates.H3&&S.gates.H4;
 else,S.gates='not applicable: no variant promoted on development data';end
end
uq_writejson(fullfile(root,'analysis',[split '_summary.json']),S);disp(jsonencode(S));
end
function T=truth(root,up,model,name)
switch model
 case {'primary','complex'},[~,ref]=load_dev(fullfile(up,'soft_refit'),model);T=ref.datasets(strcmp({ref.datasets.name},name)).target;
 case {'F1','F2'},t=load(fullfile(up,'independent',[model '_truth.mat']));T=t.tr(strcmp({t.tr.name},name)).target;
 otherwise,t=load(fullfile(root,'truth','confirm',[model '_truth.mat']));T=t.tr(strcmp({t.tr.name},name)).target;
end
end
