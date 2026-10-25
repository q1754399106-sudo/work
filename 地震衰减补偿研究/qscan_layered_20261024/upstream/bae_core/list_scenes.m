function scenes=list_scenes(root,split)
% Solver-side scene list. No truth is loaded here.
scenes=struct('model',{},'name',{},'cm',{},'D',{},'sigma',{},'referenceQ',{});
if strcmp(split,'dev')
 soft=fullfile(fileparts(root),'soft_refit_20261023');
 for which={'primary','complex'}
  [cm,ref]=load_dev(soft,which{1});
  for k=1:numel(ref.datasets)
   scenes(end+1)=struct('model',which{1},'name',ref.datasets(k).name,'cm',cm,'D',ref.datasets(k).D, ...
    'sigma',reshape(ref.datasets(k).sigma,1,[]),'referenceQ',reshape(ref.reference_Q,1,[])); %#ok<AGROW>
  end
 end
else
 if strcmp(split,'val'),val=fullfile(fileparts(root),'coherence_soft_20261023','inputs','validation');names={'V1','V2','V3'};
 else,val=fullfile(root,'inputs','independent');names={'F1','F2'};end
 for v=names
  in=load(fullfile(val,[v{1} '.mat']));assert(~isfield(in,'tr'));cm=in.config;cm.t=cm.t(:);
  for k=1:numel(in.obs)
   assert(~isfield(in.obs(k),'target'));
   scenes(end+1)=struct('model',v{1},'name',in.obs(k).name,'cm',cm,'D',in.obs(k).D, ...
    'sigma',reshape(in.obs(k).sigma,1,[]),'referenceQ',reshape(cm.Q,1,[])); %#ok<AGROW>
  end
 end
end
end
