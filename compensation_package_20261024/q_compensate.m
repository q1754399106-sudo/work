function out=q_compensate(D,cm,varargin)
% Attenuation (Q) compensation with neighbour-coherence soft shrinkage.
%
%   out = q_compensate(D, cm)                      % Standard rule, Q as given
%   out = q_compensate(D, cm, 'mode', 'qscan')     % QScan common-Q correction + Standard
%   out = q_compensate(..., 'offset', 450, 'sigma', s, 'rule', 'v1')
%
% D     nt-by-nx section gather (one common offset), all-zero columns = dead traces.
% cm    calibrated configuration: t, dt, nt, nfft, fref, fpeak, vp(1x4), h(1x3), Q(1x4).
%       The true Q is never needed; Q comes from calibration or an independent reference.
% Options
%   'mode'    'standard' (default)  FullCalRefit at the given Q, then coherence gain.
%             'qscan'               choose a common factor c on a 13-point grid in
%                                   [0.6, 2.2] by the gather-summed lasso objective,
%                                   then the same as 'standard' at Q*c.
%   'rule'    'standard' (default, lambda 0.5/16) | 'v1' (1/8) | 'soft' (lambda 1, no coherence)
%   'offset'  section offset in m (default cm.solver_offset if present, else 450)
%   'sigma'   1-by-nx noise scale; default WindowMAD on 0.02-0.25 s (supply your own if
%             that window contains signal)
%
% When to use which (see README_ZH.md): Q from a trusted independent reference ->
% 'standard'; Q only from surface-data calibration -> 'qscan'. QScan is biased towards
% under-compensation by about 10-28% on fine-layered (non-sparse) earth models, which
% costs about 3% error when the given Q is already right.
%
% out.gather       compensated section (source-wavelet domain, like FullCalRefit)
% out.unshrunk     FullCalRefit at the applied Q;  out.soft  SoftGroupRefit at the applied Q
% out.applied_Q, out.qscan_factor, out.coherence, out.lambda, out.gain, out.fallback (per trace)
% out.checks       optimisation flags, normal-equation residual, soft reproduction error
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));
% The coherence gain is sensitive to BLAS summation order (1 vs 4 threads differ by
% ~3e-3 in output); fix 4 threads, as in the archived runs, and restore afterwards.
previousThreads=maxNumCompThreads(4);restoreThreads=onCleanup(@()maxNumCompThreads(previousThreads)); %#ok<NASGU>
p=inputParser;p.addParameter('mode','standard');p.addParameter('rule','standard');p.addParameter('offset',[]);p.addParameter('sigma',[]);p.parse(varargin{:});o=p.Results;
[nt,nx]=size(D);cm.t=cm.t(:);assert(nt==numel(cm.t)&&all(isfinite(D(:))),'D must be finite, nt-by-nx with nt = numel(cm.t).');
assert(~any(isfield(cm,{'truth','true_Q'})),'Do not pass truth fields.');
offset=o.offset;if isempty(offset),if isfield(cm,'solver_offset'),offset=cm.solver_offset;else,offset=450;end,end
sigma=o.sigma;live=any(D~=0,1);
if isempty(sigma),sigma=realmin*ones(1,nx);for j=find(live),ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;end,end
sigma=reshape(sigma,1,[]);assert(numel(sigma)==nx&&all(sigma>0));
switch o.rule
 case 'standard',rule=struct('lambda_min',.5,'lambda_max',16);
 case 'v1',rule=struct('lambda_min',1,'lambda_max',8);
 case 'soft',rule=struct('lambda_min',1,'lambda_max',1);
 otherwise,error('Unknown rule %s.',o.rule);
end
factor=1;scan=[];
if strcmp(o.mode,'qscan')
 [factor,scan]=qscan_factor(D,sigma,cm,offset);
elseif ~strcmp(o.mode,'standard'),error('Unknown mode %s.',o.mode);end
used=cm;used.Q=reshape(cm.Q,1,[])*factor;d=prune_dictionary(used,offset);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;
R0=zeros(nt,nx);S0=R0;soft=R0;sup=cell(1,nx);E=sup;V=sup;flags=0;normal=0;
for b=1:11:nx
 ids=b:min(b+10,nx);dicts=repmat({d},1,numel(ids));
 ini=reflectivity_refit_prepared(D(:,ids),dicts,sigma(ids));so=reflectivity_soft_prepared(D(:,ids),dicts,sigma(ids),ini);
 R0(:,ids)=ini.reflectivity;S0(:,ids)=ini.gather;soft(:,ids)=so.gather;sup(ids)=ini.support;E(ids)=so.group_energy;V(ids)=so.group_noise_energy;
 flags=flags+ini.optimization_flags;normal=max([normal,ini.normal_residual,so.initial_normal_residual]);
end
cs=coherence_scores(R0,sup,d,live,2,.25);
s=struct('R0',R0,'supports',{sup},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'W',d.W);
one=struct('lambda_min',1,'lambda_max',1);softdiff=norm(apply_gain(s,cs.c,one)-soft,'fro')/max(norm(soft,'fro'),realmin);
assert(softdiff<1e-10,'Gain with lambda=1 must reproduce SoftGroupRefit.');
[S,~,gains]=apply_gain(s,cs.c,rule);
lam=cs.c;for j=1:nx,lam{j}=rule.lambda_min+(rule.lambda_max-rule.lambda_min)*(1-cs.c{j}).^2;lam{j}(cs.fallback{j})=1;end
out=struct('gather',S,'unshrunk',S0,'soft',soft,'applied_Q',used.Q,'qscan_factor',factor,'qscan',scan, ...
 'coherence',{cs.c},'lambda',{lam},'gain',{gains},'fallback',{cs.fallback},'sigma',sigma,'mode',o.mode,'rule',o.rule, ...
 'checks',struct('optimization_flags',flags,'max_normal_residual',normal,'soft_reproduction_error',softdiff));
end

function [factor,scan]=qscan_factor(D,sigma,used,offset)
% Selection step of bae_q_20261024/core/qscan_refit.m, unchanged.
n=size(D,1);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;grid=exp(linspace(log(.6),log(2.2),13));J=zeros(size(grid));
for q=1:numel(grid)
 c=used;c.Q=reshape(used.Q,1,[])*grid(q);d=prune_dictionary(c,offset);X=refit_lasso(d.Gram,d.B'*D(:,live),lambda(live),d.L);res=D(:,live)-d.B*X;
 J(q)=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);
end
[~,b]=min(J);factor=grid(b);scan=struct('grid',grid,'objective',J);
end
