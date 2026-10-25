function [M,S]=score_split(split)
% Truth-side scoring for 'dev', 'val' or 'ind'. Merges results/<split> and
% results/<split>_baselines. Ratios are versus FullCalRefit on the same scene.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));
scenes=list_scenes(root,split);rows={};
for k=1:numel(scenes)
 sc=scenes(k);[T,weak]=truth_of(root,split,sc);live=any(sc.D~=0,1);
 files=[];pat=sprintf('%s_Q*_%s.mat',sc.model,sc.name);
 for sub={'','_baselines','_qscan','_explore'},files=[files;dir(fullfile(root,'results',[split sub{1}],pat))];end %#ok<AGROW>
 for f=reshape(files,1,[])
  e=sscanf(regexprep(f.name,'^.*_Q([+-]\d+)_.*$','$1'),'%d');L=load(fullfile(f.folder,f.name),'g');
  for m=reshape(fieldnames(L.g),1,[])
   A=double(L.g.(m{1}));rel=norm(A(:,live)-T(:,live),'fro')/norm(T(:,live),'fro');w=NaN;
   if ~isempty(weak),w=score_weak(A,weak.reflectivity,weak.events,weak.columns,sc.cm.fpeak,sc.cm.dt);end
   rows(end+1,:)={sc.model,sc.name,e,m{1},rel,w}; %#ok<AGROW>
  end
 end
end
M=cell2table(rows,'VariableNames',{'model','condition','Q_percent','method','relative_error','weak_error'});
writetable(M,fullfile(root,'analysis',[split '_metrics_all.csv']));
methods=unique(M.method);combos=unique(M(:,{'model','condition'}));S=table();
for q=1:numel(methods)
 for c=1:height(combos)
  sel=strcmp(M.model,combos.model{c})&strcmp(M.condition,combos.condition{c});
  base=M(sel&strcmp(M.method,'FullCalRefit'),:);mine=M(sel&strcmp(M.method,methods{q}),:);
  if isempty(mine)||isempty(base),continue;end
  [~,ia,ib]=intersect(mine.Q_percent,base.Q_percent);mine=mine(ia,:);base=base(ib,:);
  wc=max(mine.relative_error);wcb=max(base.relative_error);z=mine.Q_percent==0;
  if any(z),q0=mine.relative_error(z);q0b=base.relative_error(z);w0=mine.weak_error(z)/base.weak_error(z);else,q0=NaN;q0b=NaN;w0=NaN;end
  S=[S;table(methods(q),combos.model(c),combos.condition(c),wc,wc/wcb,q0,q0/q0b,w0, ...
   'VariableNames',{'method','model','condition','wc','wc_ratio','q0','q0_ratio','weak_q0_ratio'})]; %#ok<AGROW>
 end
end
writetable(S,fullfile(root,'analysis',[split '_ratios_all.csv']));
end

function [T,weak]=truth_of(root,split,sc)
weak=[];
switch split
 case 'dev'
  soft=fullfile(fileparts(root),'soft_refit_20261023');[~,ref,input]=load_dev(soft,sc.model);
  T=ref.datasets(strcmp({ref.datasets.name},sc.name)).target;
  if strcmp(sc.model,'primary'),weak=struct('reflectivity',input.reflectivity,'events',input.event_indices_zero_based+1,'columns',[2 4 6 8]);end
 case 'val'
  t=load(fullfile(fileparts(root),'coherence_soft_20261023','truth',[sc.model '_truth.mat']));T=t.tr(strcmp({t.tr.name},sc.name)).target;
  if isfield(t,'events'),weak=struct('reflectivity',t.reflectivity,'events',t.events,'columns',t.weak);end
 case 'ind'
  t=load(fullfile(root,'truth','independent',[sc.model '_truth.mat']));T=t.tr(strcmp({t.tr.name},sc.name)).target;
end
end
