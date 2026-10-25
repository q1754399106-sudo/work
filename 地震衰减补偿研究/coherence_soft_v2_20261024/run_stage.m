function result=run_stage(stage,cands,reference,prep)
% Evaluate gain candidates on the 50 development scenes and apply the
% registered stage gate against `reference` (the current best rule).
% cands/reference: structs with fields label, opts and optional cfun(s)->c.
% prep(s,d) (optional) augments a scene once per stage (e.g. VQ); cached.
root=fileparts(mfilename('fullpath'));base=fileparts(root);v1=fullfile(base,'coherence_soft_20261023');
addpath(fullfile(root,'core'),fullfile(v1,'core'),fullfile(v1,'source'));maxNumCompThreads(4);
if nargin<4,prep=[];end
files=dir(fullfile(root,'cache','scenes','*.mat'));all=[reference,cands];rows={};started=tic;
for i=1:numel(files)
 S=load(fullfile(files(i).folder,files(i).name),'s');s=S.s;s.R0=full(s.R0);
 T=load(fullfile(root,'cache','truth',files(i).name));d=prune_dictionary(s.cm,450);s.W=d.W;
 if ~isempty(prep)
  pf=fullfile(root,'cache',['prep_' stage],files(i).name);
  if isfile(pf),P=load(pf,'s');s=P.s;s.R0=full(s.R0);s.W=d.W;
  else,s=prep(s,d);if ~isfolder(fileparts(pf)),mkdir(fileparts(pf));end;keep=rmfield(s,'W');keep.R0=sparse(keep.R0);x=s;s=keep;save(pf,'s','-v7');s=x;end
 end
 for b={'FullCalRefit','PrunedFullCalRefit','SoftGroupRefit'}
  rows{end+1}=row(s.name,b{1},score_gather(T.gathers.(b{1}),T.T,s.cm),T.T.family); %#ok<AGROW>
 end
 for c=1:numel(all)
  cc=s.c;if isfield(all(c),'cfun')&&~isempty(all(c).cfun),cc=all(c).cfun(s);end
  A=apply_gain(s,cc,all(c).opts);rows{end+1}=row(s.name,all(c).label,score_gather(A,T.T,s.cm),T.T.family); %#ok<AGROW>
 end
 if mod(i,10)==0,fprintf('[%s] %d/50 scenes %.1fs\n',stage,i,toc(started));end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis',['stage_' stage '_metrics.csv']));
ref=M(strcmp(M.method,reference.label),:);summary=struct([]);
for c=1:numel(all)
 x=M(strcmp(M.method,all(c).label),:);r=@(f)mean(x.(f)./ref.(f),'omitnan');
 q=contains(x.scene,'Q-');s=struct('label',all(c).label,'overall',r('overall'),'weak',r('weak'),'fault',r('fault'),'pinch',r('pinch'), ...
  'low_energy',r('low_energy'),'leak',r('leak'),'Q_under_overall',mean(x.overall(q)./ref.overall(q)), ...
  'four_overall',mean(x.overall(strcmp(x.family,'four'))./ref.overall(strcmp(x.family,'four'))), ...
  'complex_overall',mean(x.overall(strcmp(x.family,'complex'))./ref.overall(strcmp(x.family,'complex'))), ...
  'scenes_error_gt_1',sum(x.overall>1),'mean_overall_abs',mean(x.overall));
 s.pass=s.overall<=.99&&s.weak<=1.02&&s.fault<=1.02&&s.pinch<=1.02;
 summary=[summary,s]; %#ok<AGROW>
end
result=struct('stage',stage,'reference',reference.label,'candidates',summary);
uq_writejson(fullfile(root,'analysis',['stage_' stage '_summary.json']),result);
fprintf('%-26s %7s %7s %7s %7s %7s %7s %7s %7s %5s %s\n','label','overall','weak','fault','pinch','lowE','leak','Qunder','cplx','>1','pass');
for s=summary,fprintf('%-26s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %5d %d\n',s.label,s.overall,s.weak,s.fault,s.pinch,s.low_energy,s.leak,s.Q_under_overall,s.complex_overall,s.scenes_error_gt_1,s.pass);end
for b={'SoftGroupRefit','PrunedFullCalRefit'}
 x=M(strcmp(M.method,b{1}),:);fprintf('%-26s %7.4f %7.4f %7.4f %7.4f  (ratio to reference; >1 scenes %d)\n',b{1},mean(x.overall./ref.overall),mean(x.weak./ref.weak,'omitnan'),mean(x.fault./ref.fault,'omitnan'),mean(x.pinch./ref.pinch,'omitnan'),sum(x.overall>1));
end
end

function r=row(scene,method,m,family)
r=struct('scene',scene,'family',family,'method',method,'overall',m.overall,'weak',m.weak,'fault',m.fault,'pinch',m.pinch,'low_energy',m.low_energy,'leak',m.leak);
end
