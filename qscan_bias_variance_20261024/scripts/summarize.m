function summarize(varargin)
% Stage 3: merge replicate partials and split each estimator's error into bias^2 + variance.
% Finite-R corrections: V_u = V_ML*R/(R-1), Bias2_u = ||Xbar-T||^2 - V_u/R, hence
% Bias2_u + V_u equals the empirical MSE exactly (checked against the per-replicate scores).
P=bv_setup();nf=numel(P.factors);slots=nf+1;
labels=[arrayfun(@(c)sprintf('c=%.4f',c),P.factors,'UniformOutput',false),{'QScan'}];
stage0=readtable(fullfile(P.root,'analysis','reproduction_check.csv'),'TextType','string');
clean=readtable(fullfile(P.root,'analysis','clean_sweep.csv'),'TextType','string');
dec={};pair={};selrows={};win={};prof={};spec={};slope={};reprows={};objrows={};
which=1:numel(P.scenes);if nargin,which=[varargin{:}];end   % optional scene subset (pipeline check)
for s=which
 sc=load_scene(P,s);live=any(sc.D~=0,1);T=sc.T;ET=sum(T(:,live).^2,'all');t=sc.cm.t;
 files=dir(fullfile(P.root,'results',sprintf('rep_%s_%s_*.mat',sc.model,sc.condition)));
 M1=0;M2=0;C=0;REL=[];SUP=[];OBJ=[];MIS=[];PEN=[];ACT=[];WM=[];rep=[];seed=[];best=[];sm=[];scv=[];fl=[];nr=[];
 for f=reshape(files,1,[])
  L=load(fullfile(f.folder,f.name));M1=M1+L.M1;M2=M2+L.M2;C=C+L.C;REL=[REL;L.REL];SUP=[SUP;L.SUP]; %#ok<AGROW>
  OBJ=[OBJ;L.OBJ];MIS=[MIS;L.MIS];PEN=[PEN;L.PEN];ACT=[ACT;L.ACT];WM=[WM;L.WM];rep=[rep;L.replicate]; %#ok<AGROW>
  seed=[seed;L.seed];best=[best;L.best];sm=[sm;L.sigma_mean];scv=[scv;L.sigma_cv];fl=[fl;L.flags];nr=[nr;L.normal]; %#ok<AGROW>
  ids=L.ids;edges=L.edges;
 end
 assert(isequal(sort(rep),(1:P.replicates)'),'Replicates missing or duplicated.');R=numel(rep);
 Xbar=M1/R;V=max(M2/R-Xbar.^2,0)*R/(R-1);
 a0=stage0(stage0.model==sc.model&stage0.condition==sc.condition,:);
 cl=clean(clean.model==sc.model&clean.condition==sc.condition,:);
 cs=load(fullfile(P.root,'results',sprintf('clean_%s_%s.mat',sc.model,sc.condition)));
 mse=zeros(1,slots);b2=mse;va=mse;
 for f=1:slots
  va(f)=sum(V(:,live,f),'all');b2(f)=sum((Xbar(:,live,f)-T(:,live)).^2,'all')-va(f)/R;mse(f)=b2(f)+va(f);
  check=abs(mean(REL(:,f).^2)*ET-mse(f))/mse(f);assert(check<1e-9,'MSE identity failed.');
  if f<slots,factor=P.factors(f);cleanrel=cl.relative_error(cl.slot==f);else,factor=mode(P.grid(best));cleanrel=cl.relative_error(cl.slot==1+cs.sel.best);end
  arch=NaN;if f==1,arch=a0.rel_base;elseif f==slots,arch=a0.rel_scan;end
  dec(end+1,:)={sc.model,sc.condition,f,labels{f},factor,R,mse(f)/ET,b2(f)/ET,va(f)/ET,b2(f)/mse(f),sqrt(mse(f)/ET), ...
   mean(REL(:,f)),median(REL(:,f)),std(REL(:,f)),prctile(REL(:,f),5),prctile(REL(:,f),95),arch,cleanrel,mean(SUP(:,f)),check}; %#ok<AGROW>
 end
 ratio=REL(:,slots)./REL(:,1);dif=REL(:,slots)-REL(:,1);
 vdiff=sum(V(:,live,1)+V(:,live,slots)-2*(C(:,live)/R-Xbar(:,live,1).*Xbar(:,live,slots))*R/(R-1),'all');
 pair(end+1,:)={sc.model,sc.condition,R,mse(1)/ET,mse(slots)/ET,(mse(slots)-mse(1))/ET,(b2(slots)-b2(1))/ET,(va(slots)-va(1))/ET, ...
  (mse(slots)-mse(1))/mse(1),(b2(slots)-b2(1))/mse(1),(va(slots)-va(1))/mse(1),vdiff/ET, ...
  mean(ratio),median(ratio),prctile(ratio,5),prctile(ratio,95),mean(ratio>1),mean(dif),std(dif)/sqrt(R), ...
  a0.rel_scan/a0.rel_base,mean(best==8),numel(unique(best)),cs.sel.selected_factor}; %#ok<AGROW>
 for q=1:numel(P.grid)
  selrows(end+1,:)={sc.model,sc.condition,q,P.grid(q),sum(best==q),mean(best==q), ...
   mean(OBJ(:,q)-min(OBJ,[],2)),mean(MIS(:,q)-MIS(:,8)),mean(PEN(:,q)-PEN(:,8)),mean(ACT(:,q)), ...
   cs.sel.objective(q)-min(cs.sel.objective),cs.sel.misfit(q)-cs.sel.misfit(8),cs.sel.penalty(q)-cs.sel.penalty(8)}; %#ok<AGROW>
 end
 % Diagnostic windows (35 traces x 26 windows, same masks as theory_predict_diagnostic).
 nw=numel(edges)-1;w=0;
 for j=ids
  for k=1:nw
   w=w+1;m=t>=edges(k)-1e-9&t<edges(k+1)-1e-9;tm=T(m,j);ne=sum(m);
   bb=sum((Xbar(m,j,1)-tm).^2)-sum(V(m,j,1))/R;vb=sum(V(m,j,1));
   bs=sum((Xbar(m,j,slots)-tm).^2)-sum(V(m,j,slots))/R;vs=sum(V(m,j,slots));
   win(end+1,:)={sc.model,sc.condition,j,k,edges(k),edges(k+1),mean(tm.^2),bb/ne,vb/ne,(bb+vb)/ne,bs/ne,vs/ne,(bs+vs)/ne, ...
    mean(WM(:,w,2)>WM(:,w,1)),mean(WM(:,w,2)-WM(:,w,1)),std(WM(:,w,2)-WM(:,w,1))/sqrt(R)}; %#ok<AGROW>
  end
 end
 % Time profile over the whole gather (100 ms windows, all traces).
 pe=0:.1:1.8;
 for k=1:numel(pe)-1
  m=t>=pe(k)-1e-9&(t<pe(k+1)-1e-9|(k==numel(pe)-1&t<=pe(end)+1e-9));tw=sum(T(m,live).^2,'all');
  for f=1:slots
   vv=sum(V(m,live,f),'all');bb=sum((Xbar(m,live,f)-T(m,live)).^2,'all')-vv/R;
   cb=NaN;if f<slots,cb=sum((double(cs.X(m,live,f))-T(m,live)).^2,'all');else,cb=sum((double(cs.X(m,live,1+cs.sel.best))-T(m,live)).^2,'all');end
   prof(end+1,:)={sc.model,sc.condition,f,labels{f},pe(k),pe(k+1),tw/ET,bb/ET,vv/ET,cb/ET}; %#ok<AGROW>
  end
 end
 % Least-squares transfer from target to mean output, |H(f)|, per time band.
 bands=[.30 .80;.80 1.20;1.20 1.60];nfft=256;fr=(0:nfft/2)'/(nfft*sc.cm.dt);
 for b=1:size(bands,1)
  m=t>=bands(b,1)-1e-9&t<bands(b,2)-1e-9;hw=hann_window(sum(m));FT=fft(T(m,live).*hw,nfft);FT=FT(1:nfft/2+1,:);
  den=sum(abs(FT).^2,2);band=fr>=10&fr<=60;
  for f=1:2*slots-1
   if f<=slots,Y=Xbar(m,live,f);src='replicate_mean';lab=labels{f};else,Y=double(cs.X(m,live,f-slots));src='noise_free';lab=labels{f-slots};end
   FY=fft(Y.*hw,nfft);FY=FY(1:nfft/2+1,:);H=sum(FY.*conj(FT),2)./den;coh=abs(sum(FY.*conj(FT),2)).^2./(den.*sum(abs(FY).^2,2));
   keep=fr<=100;
   spec=[spec;[repmat({sc.model,sc.condition,src,lab,bands(b,1),bands(b,2)},sum(keep),1),num2cell([fr(keep),abs(H(keep)),coh(keep),den(keep)/max(den)])]]; %#ok<AGROW>
   wts=den(band);cf=lscov([ones(sum(band),1),fr(band)],log(abs(H(band))),wts);
   slope(end+1,:)={sc.model,sc.condition,src,lab,bands(b,1),bands(b,2),cf(1),cf(2),-cf(2)/pi}; %#ok<AGROW>
  end
 end
 for i=1:R
  reprows(end+1,:)=[{sc.model,sc.condition,rep(i),seed(i),P.grid(best(i)),sm(i),scv(i),fl(i),nr(i)},num2cell(REL(i,:)),{SUP(i,1),SUP(i,slots)}]; %#ok<AGROW>
  for q=1:numel(P.grid),objrows(end+1,:)={sc.model,sc.condition,rep(i),q,P.grid(q),OBJ(i,q),MIS(i,q),PEN(i,q),ACT(i,q)};end %#ok<AGROW>
 end
 % Rebuildable convenience copy (cache/ stays local); the partials in results/ are the record.
 Xbar=single(Xbar);V=single(V);factors=P.factors; %#ok<NASGU>
 if ~isfolder(fullfile(P.root,'cache')),mkdir(fullfile(P.root,'cache'));end
 save(fullfile(P.root,'cache',sprintf('moments_%s_%s.mat',sc.model,sc.condition)),'Xbar','V','R','factors','labels','-v7');
 fprintf('%s %s R=%d: base mse %.4f (bias2 %.4f, var %.4f); QScan mse %.4f (bias2 %.4f, var %.4f)\n',sc.model,sc.condition,R, ...
  mse(1)/ET,b2(1)/ET,va(1)/ET,mse(slots)/ET,b2(slots)/ET,va(slots)/ET);
end
A=fullfile(P.root,'analysis');
writetable(cell2table(dec,'VariableNames',{'model','condition','slot','estimator','factor','R','mse_rel','bias2_rel','var_rel', ...
 'bias_fraction','rmse_rel','mean_rel_error','median_rel_error','sd_rel_error','p05_rel_error','p95_rel_error', ...
 'archived_rel_error','noise_free_rel_error','mean_support','identity_check'}),fullfile(A,'decomposition.csv'));
writetable(cell2table(pair,'VariableNames',{'model','condition','R','mse_base','mse_qscan','delta_mse','delta_bias2','delta_var', ...
 'delta_mse_over_base','delta_bias2_over_base','delta_var_over_base','var_of_paired_difference','mean_ratio','median_ratio', ...
 'p05_ratio','p95_ratio','fraction_degraded','mean_rel_difference','se_rel_difference','archived_ratio', ...
 'fraction_selected_1p2803','distinct_selected','noise_free_selected_factor'}),fullfile(A,'qscan_vs_base.csv'));
writetable(cell2table(selrows,'VariableNames',{'model','condition','grid_index','factor','times_selected','fraction_selected', ...
 'mean_objective_above_min','mean_misfit_minus_at_1p2803','mean_penalty_minus_at_1p2803','mean_lasso_active', ...
 'noise_free_objective_above_min','noise_free_misfit_minus_at_1p2803','noise_free_penalty_minus_at_1p2803'}),fullfile(A,'qscan_selection.csv'));
writetable(cell2table(win,'VariableNames',{'model','condition','trace','window','start_s','end_s','target_ms','bias2_base','var_base', ...
 'mse_base','bias2_qscan','var_qscan','mse_qscan','p_window_degraded','mean_delta_mse','se_delta_mse'}),fullfile(A,'window_moments.csv'));
writetable(cell2table(prof,'VariableNames',{'model','condition','slot','estimator','start_s','end_s','target_energy_share', ...
 'bias2_rel','var_rel','noise_free_sq_error_rel'}),fullfile(A,'time_profile.csv'));
writetable(cell2table(spec,'VariableNames',{'model','condition','source','estimator','band_start_s','band_end_s','frequency_hz', ...
 'transfer_abs','coherence','target_power_norm'}),fullfile(A,'spectral_transfer.csv'));
writetable(cell2table(slope,'VariableNames',{'model','condition','source','estimator','band_start_s','band_end_s','log_gain_intercept', ...
 'log_gain_slope_per_hz','residual_attenuation_s'}),fullfile(A,'residual_attenuation.csv'));
relnames=strcat('rel_',strrep(strrep(labels,'=',''),'.','p'));
writetable(cell2table(reprows,'VariableNames',[{'model','condition','replicate','seed','qscan_factor','sigma_mean','sigma_cv', ...
 'optimization_flags','max_normal_residual'},relnames,{'support_base','support_qscan'}]),fullfile(A,'replicates.csv'));
writetable(cell2table(objrows,'VariableNames',{'model','condition','replicate','grid_index','factor','objective','misfit','penalty', ...
 'lasso_active_per_trace'}),fullfile(A,'replicate_qscan_objective.csv'));
fprintf('DONE stage 3\n');
end

function w=hann_window(n)
w=.5-.5*cos(2*pi*(0:n-1)'/(n-1));
end
