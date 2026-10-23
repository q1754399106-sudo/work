function run_development()
% Development round on the 20 soft_refit scenes; selects and freezes one rule.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
softroot=fullfile(fileparts(root),'soft_refit_20261023');started=tic;
plan=jsondecode(fileread(fullfile(root,'inputs','PLAN.json')));grid=plan.development.grid;
assert(~isfile(fullfile(root,'analysis','SELECTED.json')),'Development already frozen; do not rerun selection.');
test_coherence_math();
primary=jsondecode(fileread(fullfile(softroot,'inputs','primary.json')));events=primary.event_indices_zero_based+1;
rows={};lambdas=[];eventc=[];otherc=[];maxsoftdiff=0;
base=struct('half_width',2,'dip_periods',.25,'lambda_min',1,'lambda_max',1);
for protocol={'primary','complex'}
 which=protocol{1};ref=jsondecode(fileread(fullfile(softroot,'results',which,'reference.json')));
 for errorpct=[-30 -10 0 10 30]
  res=jsondecode(fileread(fullfile(softroot,'results',which,sprintf('Q_%+d.json',errorpct))));
  cm=scene_config(softroot,which,ref,errorpct);dict=prune_dictionary(cm,450);dicts=repmat({dict},1,99);
  for k=1:numel(ref.datasets)
   data=ref.datasets(k);rec=res.records(k);assert(strcmp(rec.name,data.name));D=data.D;target=data.target;
   live=any(D~=0,1);sc=load_soft_scene(softroot,which,errorpct,data.name,D);
   assert(norm(sc.soft-rec.soft,'fro')<=1e-10*norm(rec.soft,'fro'));
   same=coherence_soft_gather(sc.R0,sc.supports,dicts,sc.E,sc.V,live,base);
   softdiff=norm(same.gather-rec.soft,'fro')/norm(rec.soft,'fro');maxsoftdiff=max(maxsoftdiff,softdiff);
   assert(softdiff<1e-10,'lambda=1 does not reproduce SoftGroupRefit.');
   methods={'unpruned',rec.unpruned;'hard',rec.hard;'soft',rec.soft};
   for g=1:numel(grid)
    o=base;o.lambda_min=grid(g).lambda_min;o.lambda_max=grid(g).lambda_max;
    out=coherence_soft_gather(sc.R0,sc.supports,dicts,sc.E,sc.V,live,o);
    methods(end+1,:)={sprintf('coh_%g_%g',o.lambda_min,o.lambda_max),out.gather}; %#ok<AGROW>
    if g==1
     % lambda-independent coherence diagnostics; truth used for this diagnostic only.
     for j=find(live)
      if strcmp(which,'primary')
       cen=group_centres(sc.supports{j},dict);near=any(abs(cen-events(j,:))<=round(.5/(dict.frequency*dict.dt)),2);
       eventc=[eventc;out.coherence{j}(near)]; %#ok<AGROW>
       otherc=[otherc;out.coherence{j}(~near)]; %#ok<AGROW>
      end
     end
    end
   end
   for m=1:size(methods,1)
    A=methods{m,2};row=struct('model',which,'condition',data.name,'Q_percent',errorpct,'method',methods{m,1}, ...
     'overall_error',norm(A-target,'fro')/norm(target,'fro'),'weak_error',NaN,'weak_gain',NaN);
    if strcmp(which,'primary'),[row.weak_error,row.weak_gain]=score_weak(A,primary.reflectivity,events,[2 4 6 8],cm.fpeak,cm.dt);end
    rows{end+1}=row; %#ok<AGROW>
   end
   fprintf('%s Q%+d %s done, elapsed %.1fs\n',which,errorpct,data.name,toc(started));
  end
 end
end
T=struct2table([rows{:}]);writetable(T,fullfile(root,'analysis','development_metrics.csv'));
% Pre-registered selection.
scenes=unique(T(:,{'model','condition','Q_percent'}),'rows');summary=struct([]);
for g=1:numel(grid)
 name=sprintf('coh_%g_%g',grid(g).lambda_min,grid(g).lambda_max);ratio=zeros(height(scenes),1);weak=nan(height(scenes),1);
 for s=1:height(scenes)
  pick=@(m)T(strcmp(T.model,scenes.model{s})&strcmp(T.condition,scenes.condition{s})&T.Q_percent==scenes.Q_percent(s)&strcmp(T.method,m),:);
  a=pick(name);b=pick('soft');ratio(s)=a.overall_error/b.overall_error;
  if strcmp(scenes.model{s},'primary'),weak(s)=a.weak_error/b.weak_error;end
 end
 summary(g).name=name;summary(g).lambda_min=grid(g).lambda_min;summary(g).lambda_max=grid(g).lambda_max;
 summary(g).mean_overall_ratio=mean(ratio);summary(g).mean_weak_ratio=mean(weak,'omitnan');
 summary(g).primary_overall_ratio=mean(ratio(strcmp(scenes.model,'primary')));summary(g).complex_overall_ratio=mean(ratio(strcmp(scenes.model,'complex')));
 summary(g).feasible=summary(g).mean_weak_ratio<=1.02;summary(g).wins=sum(ratio<1);
end
feasible=find([summary.feasible]);selected=[];
if ~isempty(feasible)
 best=min([summary(feasible).mean_overall_ratio]);near=feasible([summary(feasible).mean_overall_ratio]<=best*1.005);
 key=[[summary(near).lambda_max]',-[summary(near).lambda_min]'];[~,order]=sortrows(key);selected=summary(near(order(1)));
end
uq_writejson(fullfile(root,'analysis','development_summary.json'),struct('candidates',summary,'max_lambda1_soft_difference',maxsoftdiff, ...
 'event_coherence_quantiles',quantile(eventc,[.1 .25 .5 .75 .9]),'other_coherence_quantiles',quantile(otherc,[.1 .25 .5 .75 .9]), ...
 'event_groups',numel(eventc),'other_groups',numel(otherc),'elapsed_seconds',toc(started)));
if isempty(selected),fprintf('No feasible candidate; nothing frozen.\n');return;end
% Median lambda realised by the winner on all development scenes (UniformLambda ablation).
lam=[];o=base;o.lambda_min=selected.lambda_min;o.lambda_max=selected.lambda_max;
for protocol={'primary','complex'}
 which=protocol{1};ref=jsondecode(fileread(fullfile(softroot,'results',which,'reference.json')));
 for errorpct=[-30 -10 0 10 30]
  cm=scene_config(softroot,which,ref,errorpct);dicts=repmat({prune_dictionary(cm,450)},1,99);
  for k=1:numel(ref.datasets)
   D=ref.datasets(k).D;sc=load_soft_scene(softroot,which,errorpct,ref.datasets(k).name,D);
   out=coherence_soft_gather(sc.R0,sc.supports,dicts,sc.E,sc.V,any(D~=0,1),o);
   for j=1:99,lam=[lam;out.lambda{j}(~out.fallback{j})];end %#ok<AGROW>
  end
 end
end
frozen=struct('name','CoherenceSoftRefit','half_width',2,'dip_periods',.25,'lambda_min',selected.lambda_min,'lambda_max',selected.lambda_max, ...
 'noalign_dip_periods',0,'uniform_lambda',median(lam),'selected',selected,'frozen_before_validation',true);
uq_writejson(fullfile(root,'analysis','SELECTED.json'),frozen);
fprintf('Selected %s: overall ratio %.4f, weak ratio %.4f, uniform lambda %.3f, elapsed %.1fs\n',selected.name, ...
 selected.mean_overall_ratio,selected.mean_weak_ratio,median(lam),toc(started));
end

function cm=scene_config(softroot,which,ref,errorpct)
input=jsondecode(fileread(fullfile(softroot,'inputs',[which '.json'])));
if strcmp(which,'primary')
 cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';
else
 p=input.public;for field={'vp','h','offset','lambda_grid'},f=field{1};p.(f)=reshape(p.(f),1,[]);end
 p.t=p.t(:);cm=uq_inverse_config(p,input.estimate);
end
cm.Q=reshape(ref.reference_Q,1,[])*(1+errorpct/100);
end

function cen=group_centres(support,d)
support=support(:);if isempty(support),cen=zeros(0,1);return;end
g=cumsum([1;diff(support)>max(1,round(.5/(d.frequency*d.dt)))]);cen=accumarray(g,support,[],@mean);
end
