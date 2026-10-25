function run_trueq()
% True-Q upper-bound diagnosis on V5 (inputs/PLAN_ZH.md). Truth enters only the
% oracle dictionaries and the scoring.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);started=tic;
in=load(fullfile(root,'dev','V5.mat'));TT=load(fullfile(root,'dev','V5_truth.mat'));cm=in.config;cm.t=cm.t(:);nt=cm.nt;nx=99;
model=gen_complex_truth(TT.g);
for k=1:2,dv=norm(model.datasets(k).D-in.obs(k).D,'fro')/norm(in.obs(k).D,'fro');fprintf('regenerated %s D difference %.3g\n',in.obs(k).name,dv);assert(dv<1e-10);end
cases=struct('name',{'noisy_Moderate','noisy_Strong','clean_primary'},'D',{in.obs(1).D,in.obs(2).D,model.datasets(1).primary_att}, ...
 'sigma',{in.obs(1).sigma,in.obs(2).sigma,in.obs(1).sigma},'k',{1,2,1});
for c=1:3,cases(c).D(:,~model.live)=0;cases(c).sigma=reshape(cases(c).sigma,1,[]);end
dcal=prune_dictionary(cm,450);Ws=dcal.W./reshape(dcal.scale,1,[]);dcal.source_gram=Ws'*Ws;
types={'Cal','TrueQ','TrueQW'};rows={};wavelets=zeros(2*nt,nx);
for ty=1:3
 st=repmat(struct('R0',zeros(nt,nx),'S0',zeros(nt,nx),'soft',zeros(nt,nx),'sup',{cell(1,nx)},'E',{cell(1,nx)},'V',{cell(1,nx)}),1,3);
 for b=1:9
  ids=(b-1)*11+(1:11);dicts=cell(1,11);
  for q=1:11
   j=ids(q);
   switch types{ty}
    case 'Cal',dicts{q}=dcal;
    case 'TrueQ',dicts{q}=oracle_dict(cm,model.truth.a(:,j),cm.fpeak,0,[],j);
    case 'TrueQW',[dicts{q},wavelets(:,j)]=oracle_dict(cm,model.truth.a(:,j),model.truth.fpeak(j),model.truth.phase(j),[-.16 .012;.04 .025],j);
   end
  end
  for c=1:3
   D=cases(c).D(:,ids);s=cases(c).sigma(ids);ini=reflectivity_refit_prepared(D,dicts,s);so=reflectivity_soft_prepared(D,dicts,s,ini);
   st(c).R0(:,ids)=ini.reflectivity;st(c).S0(:,ids)=ini.gather;st(c).soft(:,ids)=so.gather;st(c).sup(ids)=ini.support;st(c).E(ids)=so.group_energy;st(c).V(ids)=so.group_noise_energy;
   assert(ini.optimization_flags==0);
  end
  fprintf('%s block %d/9 %.0fs\n',types{ty},b,toc(started));
 end
 for c=1:3
  live=any(cases(c).D~=0,1);cs=coherence_scores(st(c).R0,st(c).sup,dcal,live,2,.25);
  s=struct('R0',st(c).R0,'supports',{st(c).sup},'E',{st(c).E},'V',{st(c).V},'groups',{cs.groups},'fallback',{cs.fallback},'W',dcal.W);
  [Sstd,Rstd]=apply_gain(s,cs.c,struct('lambda_min',.5,'lambda_max',16));
  if strcmp(types{ty},'TrueQW'),for j=1:nx,Sstd(:,j)=wavelet_matrix(wavelets(:,j),nt)*Rstd(:,j);end
  else,one=struct('lambda_min',1,'lambda_max',1);assert(norm(apply_gain(s,cs.c,one)-st(c).soft,'fro')<=1e-10*norm(st(c).soft,'fro'));end
  k=cases(c).k;T=struct('target',TT.tr(k).target,'masks',struct('fault',TT.regions.fault,'pinch',TT.regions.pinch,'low_energy',TT.regions.low_energy(:,:,k),'deep',TT.regions.deep));
  G=struct('FullCalRefit',st(c).S0,'SoftGroupRefit',st(c).soft,'Standard',Sstd);
  for m=fieldnames(G)'
   x=score_gather(G.(m{1}),T,cm);
   rows{end+1}=struct('dictionary',types{ty},'data',cases(c).name,'method',m{1},'overall',x.overall,'low_energy',x.low_energy,'fault',x.fault,'pinch',x.pinch,'deep',x.deep,'leak',x.leak); %#ok<AGROW>
   fprintf('%-7s %-15s %-15s overall %.4f lowE %.4f fault %.4f pinch %.4f deep %.4f\n',types{ty},cases(c).name,m{1},x.overall,x.low_energy,x.fault,x.pinch,x.deep);
  end
  if c==2||c==3,save(fullfile(root,'results',sprintf('%s_%s.mat',types{ty},cases(c).name)),'G','-v7');end
 end
end
uq_writejson(fullfile(root,'analysis','trueq_metrics.json'),[rows{:}]);
fprintf('done %.0fs\n',toc(started));
end

function [d,w]=oracle_dict(cm,a,fpeak,phase,ghost,j)
% refit_dictionary with the true per-trace attenuation a(t) and wavelet.
t=cm.t;n=numel(t);w=wavelet_vector(t,cm.dt,fpeak,phase,ghost);W=wavelet_matrix(w,n);
c=cm;G=attenuation_matrix(a,c);B=G*W;scale=sqrt(sum(B.^2,1));B=B./scale;Gram=B'*B;Gram=(Gram+Gram')/2;
% 'frequency' only sets the 0.5-period grouping gap; keep the calibrated value so that
% soft groups and coherence groups coincide (fixed after the first run failed on a
% group-count mismatch for TrueQW). W itself uses the true per-trace fpeak.
Ws=W./scale;d=struct('W',W,'B',B,'Gram',Gram,'scale',scale(:),'L',norm(Gram,1),'key',sprintf('oracle_trace_%d',j), ...
 'frequency',cm.fpeak,'dt',cm.dt,'source_gram',Ws'*Ws);
end

function w=wavelet_vector(t,dt,fpeak,phase,ghost)
% Wavelet centred at index n+1 on a 2n grid (same construction as round3 wavelet_dictionary).
n=numel(t);L=2*n;tc=(0:L-1)'*dt;z=pi*fpeak*(tc-n*dt);w=(1-2*z.^2).*exp(-z.^2);
if phase~=0
 N=2^nextpow2(4*L);F=fft(w,N);F=F(1:N/2+1);rot=exp(-1i*deg2rad(phase))*ones(size(F));rot(1)=1;rot(end)=cos(deg2rad(phase));
 Gf=F.*rot;w=real(ifft([Gf;conj(Gf(end-1:-1:2))]));w=w(1:L);
end
if ~isempty(ghost),b=w;for k=1:size(ghost,1),w=w+ghost(k,1)*interp1(tc,b,tc-ghost(k,2),'linear',0);end,end
end

function W=wavelet_matrix(w,n)
W=zeros(n,n);for i=1:n,W(:,i)=w((1:n)-i+n+1);end
end
