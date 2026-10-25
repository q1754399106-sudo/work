function run_stats(split,models)
% Solver side only: per-factor statistics on 'dev' (primary, complex, F1, F2),
% 'confirm' (new F3, F4) or 'field' (U1320A lines). Truth is never loaded here.
root=fileparts(mfilename('fullpath'));bae=fullfile(fileparts(root),'bae_q_20261024');
addpath(fullfile(root,'core'),fullfile(bae,'core'),fullfile(bae,'source'));maxNumCompThreads(4);
spec=jsondecode(fileread(fullfile(root,'SPEC.json')));grid=unique([exp(linspace(log(spec.grid_min),log(spec.grid_max),spec.grid_points)) 1]);
outdir=fullfile(root,'results',split);if ~isfolder(outdir),mkdir(outdir);end;started=tic;
if strcmp(split,'dev')
 reproduce(root,bae);scenes=[list_scenes(bae,'dev') list_scenes(bae,'ind')];errs=spec.errors_percent;
elseif strcmp(split,'confirm')
 scenes=load_dir(fullfile(root,'inputs','confirm'),{'F3','F4'});errs=spec.errors_percent;
else
 scenes=struct('model',{},'name',{},'cm',{},'D',{},'sigma',{},'referenceQ',{});
 for l={'4508-32','3677-32'}
  in=load(fullfile(fileparts(root),'field_vsp_u1320a_20261024','inputs',[l{1} '.mat']));
  scenes(end+1)=struct('model',l{1},'name','field','cm',in.cm,'D',in.D,'sigma',in.sigma,'referenceQ',in.cm.Q); %#ok<AGROW>
 end
 errs=spec.field_start_percent;
end
for k=1:numel(scenes)
 sc=scenes(k);if nargin>1&&~any(strcmp(sc.model,models)),continue;end
 offset=450;if isfield(sc.cm,'solver_offset'),offset=sc.cm.solver_offset;end
 for e=reshape(errs,1,[])
  file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',sc.model,e,sc.name));if isfile(file),continue;end
  used=sc.cm;used.Q=reshape(sc.referenceQ,1,[])*(1+e/100);st=qscan_stats(sc.D,sc.sigma,used,offset,grid);st.used_Q=used.Q;
  save(file,'st','-v7');fprintf('%s Q%+d %s flags=%d elapsed=%.1fs\n',sc.model,e,sc.name,st.optimization_flags,toc(started));
 end
end
fprintf('DONE %s\n',split);
end

function reproduce(root,bae)
% A0 on the original 13-point grid must reproduce archived QScan selections.
f=fullfile(root,'analysis','reproduction_check.json');if isfile(f),return;end
sc=list_scenes(bae,'ind');sc=sc(1);g13=exp(linspace(log(.6),log(2.2),13));res=struct('scene',{},'selected',{},'archived',{});
for e=[-50 0 50]
 used=sc.cm;used.Q=sc.referenceQ*(1+e/100);st=qscan_stats(sc.D,sc.sigma,used,0,g13);k=select_factor(st,'A0');
 A=load(fullfile(bae,'results','ind_qscan',sprintf('%s_Q%+d_%s.mat',sc.model,e,sc.name)),'info');
 res(end+1)=struct('scene',sprintf('%s_%s_Q%+d',sc.model,sc.name,e),'selected',g13(k),'archived',A.info.selected_factor); %#ok<AGROW>
end
assert(all(abs([res.selected]-[res.archived])<1e-12),'A0 does not reproduce archived QScan selections.');
uq_writejson(f,res);fprintf('Reproduction check passed.\n');
end

function scenes=load_dir(dirname,names)
scenes=struct('model',{},'name',{},'cm',{},'D',{},'sigma',{},'referenceQ',{});
for v=names
 in=load(fullfile(dirname,[v{1} '.mat']));assert(~isfield(in,'tr'));cm=in.config;cm.t=cm.t(:);
 for k=1:numel(in.obs)
  scenes(end+1)=struct('model',v{1},'name',in.obs(k).name,'cm',cm,'D',in.obs(k).D,'sigma',reshape(in.obs(k).sigma,1,[]),'referenceQ',reshape(cm.Q,1,[])); %#ok<AGROW>
 end
end
end
