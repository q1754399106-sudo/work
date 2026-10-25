function run_confirm(m,s,smoke)
% Worker (PLAN_ZH.md): model m, scene s, replicates 1..R. Truth-free decision per gather: QScan
% choice b, Z1 = lasso support at b, rule index g = b-2 (floor 1) if Z1 > threshold else b. Saves
% the Standard output (refit_state -> coherence_scores -> apply_gain, as qscan_offset_20261025/
% scripts/run_offset.m) at f=1 (B), at grid(b) (D) and at grid(max(b-2,1)) (S2). Truth is not read.
P=as_setup();sp=P.models(m);sc=P.scenes(s);mdir=fullfile(P.root,'models');
if nargin>2,P.R=1;sp.name='G1';mdir=fullfile(P.root,'cache','regen_G');end   % smoke test on a development model
M=load(fullfile(mdir,[sp.name '.mat']));cm=M.config;cm.t=cm.t(:);D0=M.clean;
Qused=reshape(cm.Q,1,[])*sc.q_scale;offset=cm.solver_offset;name=[sp.name '_' sc.tag];
dir_=fullfile(P.root,'cache','out',name);if ~isfolder(dir_),mkdir(dir_);end
ddir=fullfile(P.root,'decisions');if nargin>2,ddir=fullfile(P.root,'cache','smoke');end;if ~isfolder(ddir),mkdir(ddir);end
csv=fullfile(ddir,[name '.csv']);assert(~isfile(csv),'exists');
fid=fopen(csv,'w');fprintf(fid,'model,scene,replicate,seed,best,qscan_factor,Z1,shifted,rule_index,rule_factor,seconds\n');started=tic;
for r=1:P.R
 t0=tic;seed=P.seed(m,sc.level,r);D=D0+noise_replica(D0,sc.snr_db,seed);sig=estimate_sigmas(D,cm.t);
 sel=qscan_select(D,sig,cm,Qused,offset);b=sel.best;Z1=sel.lasso_active(b);
 shifted=Z1>P.rule.threshold;g=b;if shifted,g=max(b-P.rule.shift,1);end
 f=[1 P.grid(b) P.grid(max(b-2,1))];STD=zeros([size(D0) 3],'single');
 for k=1:3
  st=refit_state(D,sig,cm,Qused*f(k),offset);
  cs=coherence_scores(st.R0,st.supports,st.d,st.live,P.coherence.half_width,P.coherence.dip_periods);
  gs=struct('R0',st.R0,'supports',{st.supports},'E',{st.E},'V',{st.V},'groups',{cs.groups},'fallback',{cs.fallback},'W',st.d.W);
  STD(:,:,k)=single(apply_gain(gs,cs.c,P.rules.Standard));
 end
 save('-v7',fullfile(dir_,sprintf('r%03d.mat',r)),'STD','f','b','g','Z1','seed');
 fprintf(fid,'%s,%s,%d,%d,%d,%.17g,%.17g,%d,%d,%.17g,%.1f\n',sp.name,sc.tag,r,seed,b,P.grid(b),Z1,shifted,g,P.grid(g),toc(t0));
 fprintf('%s r%d: c=%.4f Z1 %.2f shift %d (%.0fs, total %.0fs)\n',name,r,P.grid(b),Z1,shifted,toc(t0),toc(started));
end
fclose(fid);fprintf('DONE %s\n',name);
end

function st=refit_state(D,sig,cm,Q,offset)
% refit_state of qscan_offset_20261025/scripts/run_offset.m (flags dropped).
c=cm;c.Q=Q;d=prune_dictionary(c,offset);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);nx=size(D,2);
st=struct('R0',zeros(size(D)),'S0',zeros(size(D)),'soft',zeros(size(D)),'supports',{cell(1,nx)},'E',{cell(1,nx)},'V',{cell(1,nx)},'live',any(D~=0,1));
for b=1:ceil(nx/11)
 ids=(b-1)*11+1:min(b*11,nx);ini=reflectivity_refit_prepared(D(:,ids),dicts(1:numel(ids)),sig(ids));so=reflectivity_soft_prepared(D(:,ids),dicts(1:numel(ids)),sig(ids),ini);
 st.R0(:,ids)=ini.reflectivity;st.S0(:,ids)=ini.gather;st.soft(:,ids)=so.gather;st.supports(ids)=ini.support;st.E(ids)=so.group_energy;st.V(ids)=so.group_noise_energy;
end
st.d=d;
end
