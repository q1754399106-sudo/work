function run_combo(sets)
% QScan + coherence soft shrinkage (inputs/PLAN_QSCAN_COH.json). Solver side
% loads only dev/<set>.mat; truth is opened by score_combo.m.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
if nargin<1,sets={'V4','F1','F2','V5'};end
rules=struct('v1',struct('lambda_min',1,'lambda_max',8),'Standard',struct('lambda_min',.5,'lambda_max',16));
cdir=fullfile(root,'cache');rdir=fullfile(root,'results');if ~isfolder(rdir),mkdir(rdir);end
started=tic;flags=0;maxnormal=0;
for v=sets
 in=load(fullfile(root,'dev',[v{1} '.mat']));assert(~isfield(in,'tr')&&~isfield(in.obs,'target'));
 cm=in.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);offset=450;if isfield(cm,'solver_offset'),offset=cm.solver_offset;end
 for e=[-50 -30 -10 0 10 30 50]
  for k=1:numel(in.obs)
   out=fullfile(rdir,sprintf('%s_Q%+d_%s.mat',v{1},e,in.obs(k).name));if isfile(out),continue;end
   D=in.obs(k).D;sig=reshape(in.obs(k).sigma,1,[]);used=cm;used.Q=Q0*(1+e/100);
   [factor,J,grid]=qscan_factor(D,sig,used,offset);scanned=used;scanned.Q=used.Q*factor;
   gathers=struct();scan=struct('factor',factor,'objective',J,'grid',grid); %#ok<NASGU>
   for route={'nominal','qscan'}
    if strcmp(route{1},'nominal'),c=used;else,c=scanned;end
    [st,f,nr]=refit_state(D,sig,c,offset,fullfile(cdir,sprintf('%s_Q%+d_%s_%s',v{1},e,in.obs(k).name,route{1})));
    flags=flags+f;maxnormal=max(maxnormal,nr);
    cs=coherence_scores(st.R0,st.supports,st.d,st.live,2,.25);
    s=struct('R0',st.R0,'supports',{st.supports},'E',{st.E},'V',{st.V},'groups',{cs.groups},'fallback',{cs.fallback},'W',st.d.W);
    one=struct('lambda_min',1,'lambda_max',1);assert(norm(apply_gain(s,cs.c,one)-st.soft,'fro')<=1e-10*max(norm(st.soft,'fro'),realmin));
    if strcmp(route{1},'nominal')
     gathers.FullCalRefit=st.S0;gathers.SoftGroupRefit=st.soft;gathers.Standard=apply_gain(s,cs.c,rules.Standard);
    else
     gathers.QScanRefit=st.S0;gathers.QScan_Soft=st.soft;gathers.QScan_v1=apply_gain(s,cs.c,rules.v1);gathers.QScan_Standard=apply_gain(s,cs.c,rules.Standard);
    end
   end
   save(out,'gathers','scan','used','scanned','-v7');
   fprintf('%s Q%+d %s factor %.3f done %.1fs\n',v{1},e,in.obs(k).name,factor,toc(started));
  end
 end
end
uq_writejson(fullfile(root,'analysis',sprintf('combo_checks_%s.json',strjoin(sets,'_'))),struct('optimization_flags',flags,'max_normal_residual',maxnormal,'elapsed_seconds',toc(started)));
end

function [factor,J,grid]=qscan_factor(D,sigma,used,offset)
% Selection step of bae_q_20261024/core/qscan_refit.m, unchanged.
n=size(D,1);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;grid=exp(linspace(log(.6),log(2.2),13));J=zeros(size(grid));
for q=1:numel(grid)
 c=used;c.Q=used.Q*grid(q);d=prune_dictionary(c,offset);X=refit_lasso(d.Gram,d.B'*D(:,live),lambda(live),d.L);res=D(:,live)-d.B*X;
 J(q)=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);
end
[~,b]=min(J);factor=grid(b);
end

function [st,flags,normal]=refit_state(D,sig,c,offset,cachefile)
d=prune_dictionary(c,offset);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);nx=size(D,2);
st=struct('R0',zeros(size(D)),'S0',zeros(size(D)),'soft',zeros(size(D)),'supports',{cell(1,nx)},'E',{cell(1,nx)},'V',{cell(1,nx)},'live',any(D~=0,1));
if isfile([cachefile '.mat']),x=load([cachefile '.mat'],'st','flags','normal');st=x.st;flags=x.flags;normal=x.normal;st.d=d;return;end
flags=0;normal=0;
for b=1:ceil(nx/11)
 ids=(b-1)*11+1:min(b*11,nx);ini=reflectivity_refit_prepared(D(:,ids),dicts(1:numel(ids)),sig(ids));so=reflectivity_soft_prepared(D(:,ids),dicts(1:numel(ids)),sig(ids),ini);
 st.R0(:,ids)=ini.reflectivity;st.S0(:,ids)=ini.gather;st.soft(:,ids)=so.gather;st.supports(ids)=ini.support;st.E(ids)=so.group_energy;st.V(ids)=so.group_noise_energy;
 flags=flags+ini.optimization_flags;normal=max([normal,ini.normal_residual,so.initial_normal_residual]);
end
save([cachefile '.mat'],'st','flags','normal','-v7');st.d=d;
end
