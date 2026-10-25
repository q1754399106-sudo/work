function run_layer(split,models)
% Solver side of the depth-segmented study, entirely from the upstream/ snapshot.
% split: 'dev' (primary, complex, F1, F2), 'confirm' (new F5, F6) or 'field'.
root=fileparts(mfilename('fullpath'));up=fullfile(root,'upstream');
addpath(fullfile(root,'core'),fullfile(up,'bae_core'),fullfile(up,'bae_source'));maxNumCompThreads(4);
spec=jsondecode(fileread(fullfile(root,'SPEC.json')));outdir=fullfile(root,'results',split);if ~isfolder(outdir),mkdir(outdir);end
[scenes,errs]=scene_list(root,up,split,spec);started=tic;
for k=1:numel(scenes)
 sc=scenes(k);if nargin>1&&~any(strcmp(sc.model,models)),continue;end
 for e=reshape(errs,1,[])
  file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',sc.model,e,sc.name));if isfile(file),continue;end
  used=sc.cm;used.Q=reshape(sc.referenceQ,1,[])*(1+e/100);
  if used.Q(1)>1e5,used.Q(1)=sc.referenceQ(1);end            % field: keep the water layer's Q
  r=layer_scan(sc.D,sc.sigma,used,sc.offset,spec);save(file,'r','-v7');
  fprintf('%s Q%+d %s c0=%.3f %s evals=%d elapsed=%.1fs\n',sc.model,e,sc.name,r.c0, ...
   strjoin(arrayfun(@(v)sprintf('%s=[%s]',v.name,num2str(v.factors,'%.2f ')),r.variants,'UniformOutput',false),' '),r.evaluations,toc(started));
 end
end
fprintf('DONE %s\n',split);
end

function [scenes,errs]=scene_list(root,up,split,spec)
scenes=struct('model',{},'name',{},'cm',{},'D',{},'sigma',{},'referenceQ',{},'offset',{});errs=spec.errors_percent;
switch split
 case 'dev'
  for w={'primary','complex'}
   [cm,ref]=load_dev(fullfile(up,'soft_refit'),w{1});
   for k=1:numel(ref.datasets),d=ref.datasets(k);scenes(end+1)=mk(w{1},d.name,cm,d.D,d.sigma,cm.Q,450);end %#ok<AGROW>
  end
  scenes=[scenes from_dir(fullfile(up,'independent'),{'F1','F2'})];
 case 'confirm',scenes=from_dir(fullfile(root,'inputs','confirm'),{'F5','F6'});
 case 'field'
  for l={'4508-32','3677-32'},in=load(fullfile(up,'field',[l{1} '.mat']));scenes(end+1)=mk(l{1},'field',in.cm,in.D,in.sigma,in.cm.Q,0);end %#ok<AGROW>
  errs=spec.field_start_percent;
end
end
function s=from_dir(d,names)
s=struct('model',{},'name',{},'cm',{},'D',{},'sigma',{},'referenceQ',{},'offset',{});
for v=names
 in=load(fullfile(d,[v{1} '.mat']));assert(~isfield(in,'tr'));cm=in.config;cm.t=cm.t(:);
 for k=1:numel(in.obs),s(end+1)=mk(v{1},in.obs(k).name,cm,in.obs(k).D,in.obs(k).sigma,cm.Q,cm.solver_offset);end %#ok<AGROW>
end
end
function s=mk(model,name,cm,D,sigma,Q,offset)
s=struct('model',model,'name',name,'cm',cm,'D',D,'sigma',reshape(sigma,1,[]),'referenceQ',reshape(Q,1,[]),'offset',offset);
end
