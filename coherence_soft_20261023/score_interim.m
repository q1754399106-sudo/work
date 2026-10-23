function score_interim(root)
% Interim copy of score_validation: skips scenes not yet solved; not an acceptance result.
if nargin<1,root=fileparts(mfilename('fullpath'));end
addpath(fullfile(root,'core'));
methods={'FullCalRefit','PrunedFullCalRefit','SoftGroupRefit','CoherenceSoftRefit','NoAlign','UniformLambda'};rows={};
for name={'V1','V2','V3'}
 v=name{1};T=load(fullfile(root,'truth',[v '_truth.mat']));in=load(fullfile(root,'inputs','validation',[v '.mat']),'config');cm=in.config;
 for errorpct=[-30 -10 0 10 30]
  for k=1:numel(T.tr)
   rf=fullfile(root,'results',sprintf('%s_Q%+d_%s.mat',v,errorpct,T.tr(k).name));if ~isfile(rf),continue;end
   r=load(rf,'gathers');tg=T.tr(k).target;
   L=sqrt(movmean(tg.^2,21,1));quiet=L<.02*max(L(:));
   for m=1:numel(methods)
    A=r.gathers.(methods{m});row=struct('set',v,'condition',T.tr(k).name,'Q_percent',errorpct,'method',methods{m}, ...
     'overall_error',rel(A,tg,true(size(tg))),'background_leakage',norm(A(quiet))/norm(tg(:)), ...
     'weak_error',NaN,'weak_fault_error',NaN,'weak_pinch_error',NaN,'fault_error',NaN,'pinch_error',NaN,'low_energy_error',NaN,'deep_error',NaN);
    if isfield(T,'events')
     row.weak_error=score_weak(A,T.reflectivity,T.events,T.weak,cm.fpeak,cm.dt);
     if any(T.regions.fault(:))
      row.weak_fault_error=score_weak(A,T.reflectivity,T.events,T.weak,cm.fpeak,cm.dt,46:54);
      row.weak_pinch_error=score_weak(A,T.reflectivity,T.events,6,cm.fpeak,cm.dt,60:80);
     end
    end
    if any(T.regions.fault(:)),row.fault_error=rel(A,tg,T.regions.fault);row.pinch_error=rel(A,tg,T.regions.pinch);end
    if isfield(T.regions,'low_energy'),row.low_energy_error=rel(A,tg,T.regions.low_energy(:,:,k));row.deep_error=rel(A,tg,T.regions.deep);end
    rows{end+1}=row; %#ok<AGROW>
   end
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis','interim_metrics.csv'));
cand='CoherenceSoftRefit';base='SoftGroupRefit';ratios=struct();
for metric={'overall_error','weak_error','weak_fault_error','weak_pinch_error','fault_error','pinch_error','low_energy_error','deep_error','background_leakage'}
 f=metric{1};for other={'SoftGroupRefit','PrunedFullCalRefit','FullCalRefit','NoAlign','UniformLambda'}
  ratios.(f).(other{1})=paired(M,cand,other{1},f);
 end
end
ov=paired(M,cand,base,'overall_error');wk=paired(M,cand,base,'weak_error');ft=paired(M,cand,base,'fault_error');pc=paired(M,cand,base,'pinch_error');
under=M(strcmp(M.method,cand)&M.Q_percent<0,:);ub=M(strcmp(M.method,base)&M.Q_percent<0,:);
qunder=mean(under.overall_error./ub.overall_error);
gates=struct('A1_overall',struct('value',ov.mean_ratio,'threshold',.95,'pass',ov.mean_ratio<=.95), ...
 'A2_weak',struct('value',wk.mean_ratio,'threshold',1.02,'pass',wk.mean_ratio<=1.02), ...
 'A3_fault',struct('value',ft.mean_ratio,'threshold',1.02,'pass',ft.mean_ratio<=1.02), ...
 'A4_pinchout',struct('value',pc.mean_ratio,'threshold',1.02,'pass',pc.mean_ratio<=1.02), ...
 'A5_Q_under',struct('value',qunder,'threshold',1,'pass',qunder<=1));
names=fieldnames(gates);allpass=all(cellfun(@(n)gates.(n).pass,names));
bySet=struct();
for name={'V1','V2','V3'}
 s=M(strcmp(M.set,name{1}),:);bySet.(name{1})=struct('overall_vs_soft',paired(s,cand,base,'overall_error'),'weak_vs_soft',paired(s,cand,base,'weak_error'), ...
  'overall_vs_hard',paired(s,cand,'PrunedFullCalRefit','overall_error'),'weak_vs_hard',paired(s,cand,'PrunedFullCalRefit','weak_error'));
end
uq_writejson(fullfile(root,'analysis','interim_summary.json'),struct('gates',gates,'all_pass',allpass,'ratios',ratios,'by_set',bySet,'scenes',height(M)/numel(methods)));
fprintf('A1 %.4f A2 %.4f A3 %.4f A4 %.4f A5 %.4f all_pass=%d\n',ov.mean_ratio,wk.mean_ratio,ft.mean_ratio,pc.mean_ratio,qunder,allpass);
end

function e=rel(A,T,mask)
e=norm(A(mask)-T(mask))/norm(T(mask));
end

function s=paired(M,a,b,f)
x=M(strcmp(M.method,a),:);y=M(strcmp(M.method,b),:);r=x.(f)./y.(f);ok=isfinite(r);
s=struct('mean_ratio',mean(r(ok)),'wins',sum(r(ok)<1),'n',sum(ok),'min_ratio',min([r(ok);NaN]),'max_ratio',max([r(ok);NaN]));
if ~any(ok),s.mean_ratio=NaN;end
end
