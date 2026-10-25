function run_fixed(sets)
% Earliest FixedSplit (fixed_v1_*) on round3 data, same Q-error grid and solver
% offset as run_combo.m. Solver side loads only dev/<set>.mat; truth is opened
% by score_fixed.m.
root=fileparts(mfilename('fullpath'));r3=fullfile(fileparts(root),'round3_20261024');
addpath(fullfile(r3,'core'),fullfile(r3,'source'));maxNumCompThreads(4);
if nargin<1,sets={'F3','F4','V5'};end
rdir=fullfile(root,'results');if ~isfolder(rdir),mkdir(rdir);end
started=tic;maxnormal=0;maxreplay=0;
for v=sets
 in=load(fullfile(r3,'dev',[v{1} '.mat']));assert(~isfield(in,'tr')&&~isfield(in.obs,'target'));
 cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);offset=450;if isfield(cm,'solver_offset'),offset=cm.solver_offset;end
 for e=[-50 -30 -10 0 10 30 50]
  used=cm;used.Q=Q0*(1+e/100);op=[];
  for k=1:numel(in.obs)
   out=fullfile(rdir,sprintf('%s_Q%+d_%s.mat',v{1},e,in.obs(k).name));if isfile(out),continue;end
   D=in.obs(k).D;nx=size(D,2);
   if isempty(op),op=fixed_v1_operators(used,repmat(offset,1,nx));end
   [S,normal,replay]=fixed_v1_apply(D,op);maxnormal=max(maxnormal,normal);maxreplay=max(maxreplay,replay);
   save(out,'S','used','offset','normal','replay','-v7');
   fprintf('%s Q%+d %s done normal %.2g %.1fs\n',v{1},e,in.obs(k).name,normal,toc(started));
  end
 end
end
uq_writejson(fullfile(root,'analysis',sprintf('fixed_checks_%s.json',strjoin(sets,'_'))), ...
 struct('max_normal_residual',maxnormal,'max_replay_difference',maxreplay,'elapsed_seconds',toc(started)));
end
