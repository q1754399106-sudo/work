function prepare_dev()
% Builds the 50-scene development cache from existing checkpoints (no new solves).
% Solver-side state -> cache/scenes, truth-side scoring info -> cache/truth.
root=fileparts(mfilename('fullpath'));base=fileparts(root);
v1=fullfile(base,'coherence_soft_20261023');softroot=fullfile(base,'soft_refit_20261023');
addpath(fullfile(root,'core'),fullfile(v1,'core'),fullfile(v1,'source'));maxNumCompThreads(4);
sdir=fullfile(root,'cache','scenes');tdir=fullfile(root,'cache','truth');
if ~isfolder(sdir),mkdir(sdir);end;if ~isfolder(tdir),mkdir(tdir);end
primary=jsondecode(fileread(fullfile(softroot,'inputs','primary.json')));started=tic;maxdiff=0;
for which={'primary','complex'}
 w=which{1};ref=jsondecode(fileread(fullfile(softroot,'results',w,'reference.json')));
 input=jsondecode(fileread(fullfile(softroot,'inputs',[w '.json'])));
 if strcmp(w,'primary'),cm=input.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';
 else,p=input.public;for f={'vp','h','offset','lambda_grid'},p.(f{1})=reshape(p.(f{1}),1,[]);end;p.t=p.t(:);cm=uq_inverse_config(p,input.estimate);end
 for e=[-30 -10 0 10 30]
  res=jsondecode(fileread(fullfile(softroot,'results',w,sprintf('Q_%+d.json',e))));
  used=cm;used.Q=reshape(ref.reference_Q,1,[])*(1+e/100);
  for k=1:2
   data=ref.datasets(k);sc=load_soft_scene(softroot,w,e,data.name,data.D);
   name=sprintf('S_%s_Q%+d_%s',w,e,data.name);
   base_g=struct('FullCalRefit',res.records(k).unpruned,'PrunedFullCalRefit',res.records(k).hard,'SoftGroupRefit',res.records(k).soft);
   T=struct('target',data.target,'family',ternary(strcmp(w,'primary'),'four','complex'),'masks',struct());
   if strcmp(w,'primary'),T.reflectivity=primary.reflectivity;T.events=primary.event_indices_zero_based+1;T.weak=[2 4 6 8];end
   maxdiff=max(maxdiff,save_scene(sdir,tdir,name,used,data.D,reshape(data.sigma,1,[]),sc,base_g,T));
   fprintf('%s %.1fs\n',name,toc(started));
  end
 end
end
for v={'V1','V2','V3'}
 in=load(fullfile(v1,'inputs','validation',[v{1} '.mat']));TT=load(fullfile(v1,'truth',[v{1} '_truth.mat']));
 for e=[-30 -10 0 10 30]
  used=in.config;used.t=used.t(:);used.Q=reshape(in.config.Q,1,[])*(1+e/100);
  for k=1:2
   D=in.obs(k).D;nm=in.obs(k).name;R0=zeros(size(D));sup=cell(1,99);E=sup;V=sup;
   for b=1:9
    ids=(b-1)*11+(1:11);c=load(fullfile(v1,'cache',sprintf('%s_Q%+d_%s_%02d.mat',v{1},e,nm,b)),'initial','softout');
    R0(:,ids)=c.initial.reflectivity;sup(ids)=c.initial.support;E(ids)=c.softout.group_energy;V(ids)=c.softout.group_noise_energy;
   end
   r=load(fullfile(v1,'results',sprintf('%s_Q%+d_%s.mat',v{1},e,nm)),'gathers');
   sc=struct('R0',R0,'supports',{sup},'E',{E},'V',{V});
   T=struct('target',TT.tr(k).target,'family',ternary(strcmp(v{1},'V3'),'complex','four'),'masks',struct());
   if isfield(TT,'events'),T.reflectivity=TT.reflectivity;T.events=TT.events;T.weak=TT.weak;end
   if any(TT.regions.fault(:)),T.masks.fault=TT.regions.fault;T.masks.pinch=TT.regions.pinch;end
   if isfield(TT.regions,'low_energy'),T.masks.low_energy=TT.regions.low_energy(:,:,k);end
   name=sprintf('%s_Q%+d_%s',v{1},e,nm);
   maxdiff=max(maxdiff,save_scene(sdir,tdir,name,used,D,reshape(in.obs(k).sigma,1,[]),sc,r.gathers,T));
   fprintf('%s %.1fs\n',name,toc(started));
  end
 end
end
fprintf('Prepared 50 scenes; max |v1 recomputed - v1 stored| = %.3g (V scenes)\n',maxdiff);
end

function diffv=save_scene(sdir,tdir,name,cm,D,sigma,sc,gathers,T)
d=prune_dictionary(cm,450);live=any(D~=0,1);
cs=coherence_scores(sc.R0,sc.supports,d,live,2,.25); %#ok<NASGU>
s=struct('R0',sparse(sc.R0),'supports',{sc.supports},'E',{sc.E},'V',{sc.V},'groups',{cs.groups},'fallback',{cs.fallback}, ...
 'c',{cs.c},'lags',{cs.lags},'scores',{cs.scores},'centers',{cs.centers},'cm',cm,'D',D,'sigma',sigma,'live',live,'name',name); %#ok<NASGU>
v1=struct('lambda_min',1,'lambda_max',8);s.W=d.W;s.R0=full(s.R0);[S1,~,~]=apply_gain(s,cs.c,v1);
diffv=0;if isfield(gathers,'CoherenceSoftRefit'),diffv=norm(S1-gathers.CoherenceSoftRefit,'fro')/norm(gathers.CoherenceSoftRefit,'fro');assert(diffv<1e-10);end
gathers.CoherenceV1=S1;s=rmfield(s,'W');s.R0=sparse(s.R0);
save(fullfile(sdir,[name '.mat']),'s','-v7');save(fullfile(tdir,[name '.mat']),'T','gathers','-v7');
end

function r=ternary(c,a,b)
if c,r=a;else,r=b;end
end
