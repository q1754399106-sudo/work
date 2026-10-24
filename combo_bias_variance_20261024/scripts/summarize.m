function summarize(varargin)
% Stage 1 scoring (DESIGN_ZH.md): MSE = bias^2 + variance for every estimator, region and time
% window; the 2x2 factorial effects of QScan without (raw = C-A) and after coherence
% (coh = D-B) and their interaction; amplitude gain on strong reflections; batch standard
% errors (5 batches); J1-J5. Optional argument: scene subset (pipeline check on partial data).
P=cb_setup();A=fullfile(P.root,'analysis');which=1:numel(P.scenes);if nargin,which=[varargin{:}];end
M=P.methods;nm=numel(M);nb=5;dec={};fac={};amp={};reps={};
for s=which
 sc=load_scene(P,s);T=sc.T;t=sc.cm.t;ET=sum(T(:).^2);[n,nx]=size(T);
 files=arrayfun(@(r)fullfile(P.root,'cache','outputs',sc.name,sprintf('r%03d.mat',r)),1:sc.R,'UniformOutput',false);
 have=cellfun(@isfile,files);if ~all(have),fprintf('%s: %d/%d replicates, skipped\n',sc.name,sum(have),sc.R);continue;end
 R=sc.R;X=zeros(n,nx,R,nm,'single');factor=zeros(R,1);
 for r=1:R
  L=load(files{r},'G','info','seed');assert(L.seed==sc.seed_base+r);factor(r)=L.info.qscan_factor;
  for m=1:nm,X(:,:,r,m)=L.G.(M{m});end
 end
 % Regions from the target only (score_gather envelope) and time windows; masks partition the gather.
 Lt=sqrt(movmean(T.^2,21,1));mx=max(Lt(:));
 regions={'all',true(n,nx);'strong',Lt>=.25*mx;'weak',Lt>=.02*mx&Lt<.25*mx;'quiet',Lt<.02*mx};
 edges=[0 .30 .80 1.20 Inf];wn={'t0.00-0.30','t0.30-0.80','t0.80-1.20','t1.20-end'};
 for w=1:4,regions(end+1,:)={wn{w},repmat(t>=edges(w)-1e-9&t<edges(w+1)-1e-9,1,nx)};end %#ok<AGROW>
 strong=regions{2,2};batch=reshape(1:R,[],nb);
 % est(m,g,k,:) = [mse bias2 var] / ET for method m, region g, sample k (k=1 all replicates, 1+b batch b).
 est=zeros(nm,size(regions,1),1+nb,3);
 for m=1:nm
  Y=double(X(:,:,:,m));
  for k=1:1+nb
   if k==1,ids=1:R;else,ids=batch(:,k-1)';end
   Rk=numel(ids);Yb=mean(Y(:,:,ids),3);V=var(Y(:,:,ids),0,3);
   for g=1:size(regions,1)
    G=regions{g,2};va=sum(V(G));b2=sum((Yb(G)-T(G)).^2)-va/Rk;est(m,g,k,:)=[b2+va b2 va]/ET;
   end
  end
  % Identity check against per-replicate squared errors (all replicates, whole gather).
  sq=squeeze(sum((Y-T).^2,[1 2]))/ET;assert(abs(mean(sq)-est(m,1,1,1))<=1e-9*est(m,1,1,1),'MSE identity failed.');
  gr=squeeze(sum(Y.*T.*strong,[1 2]))/sum(T(strong).^2);gm=sum(mean(Y,3).*T.*strong,'all')/sum(T(strong).^2);
  amp(end+1,:)={sc.name,M{m},gm,mean(gr),std(gr)}; %#ok<AGROW>
  for r=1:R,reps(end+1,:)={sc.name,r,sc.seed_base+r,factor(r),M{m},sqrt(sq(r))};end %#ok<AGROW>
  for g=1:size(regions,1)
   se=std(squeeze(est(m,g,2:end,:)),0,1)/sqrt(nb);
   dec(end+1,:)=[{sc.name,sc.model,sc.snr_db,round(100*sc.q_scale)-100,R,M{m},regions{g,1},sum(T(regions{g,2}).^2)/ET}, ...
    num2cell(squeeze(est(m,g,1,:))'),num2cell(se)]; %#ok<AGROW>
  end
 end
 % Factorial effects: raw = C - A, coh = D - B, interaction = coh - raw.
 ix=@(name)find(strcmp(M,name));a=ix('FullCalRefit');b=ix('Standard');c=ix('QScanRefit');d=ix('QScan_Standard');
 for g=1:size(regions,1)
  raw=squeeze(est(c,g,:,:)-est(a,g,:,:));coh=squeeze(est(d,g,:,:)-est(b,g,:,:));int=coh-raw;
  fac(end+1,:)=[{sc.name,sc.model,sc.snr_db,round(100*sc.q_scale)-100,R,regions{g,1}},num2cell(raw(1,:)),num2cell(coh(1,:)),num2cell(int(1,:)), ...
   num2cell(std(raw(2:end,:),0,1)/sqrt(nb)),num2cell(std(coh(2:end,:),0,1)/sqrt(nb)),num2cell(std(int(2:end,:),0,1)/sqrt(nb)), ...
   {mean(factor),mode(factor),numel(unique(factor))}]; %#ok<AGROW>
 end
 fprintf('%s R=%d: A %.4f B %.4f C %.4f D %.4f | dMSE raw %+.4f coh %+.4f (dVar raw %+.4f coh %+.4f)\n',sc.name,R, ...
  est(a,1,1,1),est(b,1,1,1),est(c,1,1,1),est(d,1,1,1),est(c,1,1,1)-est(a,1,1,1),est(d,1,1,1)-est(b,1,1,1),est(c,1,1,3)-est(a,1,1,3),est(d,1,1,3)-est(b,1,1,3));
end
if isempty(dec),fprintf('No complete scene.\n');return;end
dv={'scene','model','snr_db','Q_percent','R','method','region','target_energy_share','mse','bias2','var','se_mse','se_bias2','se_var'};
writetable(cell2table(dec,'VariableNames',dv),fullfile(A,'decomposition.csv'));
fv={'scene','model','snr_db','Q_percent','R','region','raw_dmse','raw_dbias2','raw_dvar','coh_dmse','coh_dbias2','coh_dvar', ...
 'int_dmse','int_dbias2','int_dvar','se_raw_dmse','se_raw_dbias2','se_raw_dvar','se_coh_dmse','se_coh_dbias2','se_coh_dvar', ...
 'se_int_dmse','se_int_dbias2','se_int_dvar','mean_factor','mode_factor','distinct_factors'};
F=cell2table(fac,'VariableNames',fv);writetable(F,fullfile(A,'factorial.csv'));
G=cell2table(amp,'VariableNames',{'scene','method','gain_of_mean','mean_replicate_gain','sd_replicate_gain'});writetable(G,fullfile(A,'amplitude_strong.csv'));
writetable(cell2table(reps,'VariableNames',{'scene','replicate','seed','qscan_factor','method','rel_error'}),fullfile(A,'replicates.csv'));
decide(P,F,G,A);
end

function decide(P,F,G,A)
% J1-J5 exactly as written in the frozen plan; a judgement needs all of its scenes.
name=@(ids)arrayfun(@(i)P.scenes(i).name,ids,'UniformOutput',false);all_=F(strcmp(F.region,'all'),:);
pick=@(ids)all_(ismember(all_.scene,name(ids)),:);out=struct();
q0=pick(1:4);
if height(q0)==4
 rhoV=q0.coh_dvar./q0.raw_dvar;rhoB=q0.coh_dbias2./q0.raw_dbias2;
 out.J1=struct('rho_V',rhoV','verdict',ternary(sum(rhoV<.5)>=3,'substitute',ternary(sum(rhoV>=.5)>=3,'complementary','mixed')));
 out.J2=struct('rho_B',rhoB','verdict',ternary(sum(rhoB>=.5)>=3,'bias cost retained','not retained'));
 sh=struct();Fq=F(ismember(F.scene,name(1:4)),:);tot=sum(Fq.coh_dmse(strcmp(Fq.region,'all')));
 for r={'strong','weak','quiet','t0.00-0.30','t0.30-0.80','t0.80-1.20','t1.20-end'}
  sh.(matlab.lang.makeValidName(r{1}))=sum(Fq.coh_dmse(strcmp(Fq.region,r{1})))/tot;
 end
 gB=G.gain_of_mean(ismember(G.scene,name(1:4))&strcmp(G.method,'Standard'));gD=G.gain_of_mean(ismember(G.scene,name(1:4))&strcmp(G.method,'QScan_Standard'));
 out.J5=struct('share_of_coh_dmse',sh,'gain_B',gB','gain_D',gD','verdict',ternary(sh.strong>.5&&mean(gD)<mean(gB),'strong-reflection amplitude under-compensated','not mainly strong-reflection amplitude'));
end
qe=pick(5:8);
if height(qe)==4,out.J3=struct('n_pass',sum(qe.coh_dmse<0&qe.coh_dbias2<0),'verdict',ternary(sum(qe.coh_dmse<0&qe.coh_dbias2<0)>=3,'works by reducing bias','not established'));end
hs=pick(9:10);
if height(hs)==2,out.J4=struct('coh_dmse',hs.coh_dmse','verdict',ternary(all(hs.coh_dmse>0),'QScan is pure cost at high SNR','not established'));end
fid=fopen(fullfile(A,'decisions.json'),'w');fprintf(fid,'%s',jsonencode(out,'PrettyPrint',true));fclose(fid);
disp(jsonencode(out,'PrettyPrint',true));
end

function v=ternary(c,a,b)
if c,v=a;else,v=b;end
end
