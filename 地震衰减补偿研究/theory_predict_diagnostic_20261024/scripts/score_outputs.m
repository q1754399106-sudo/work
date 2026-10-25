function score_outputs()
% Truth is first made available here, after predictions are externally sealed.
root=fileparts(fileparts(mfilename('fullpath')));repo=fileparts(root);
addpath(fullfile(repo,'bae_q_20261024','core'));
assert(isfile(fullfile(root,'predictions','SEALED_SHA256.json')));
P=readtable(fullfile(root,'predictions','risk.csv'),'TextType','string');
P.error_rms=zeros(height(P),1);P.base_mse=P.error_rms;P.scan_mse=P.error_rms;P.target_rms=P.error_rms;
files=dir(fullfile(root,'inputs','*.mat'));summary={};
for fi=1:numel(files)
 s=load(fullfile(files(fi).folder,files(fi).name));weak=[];
 if any(strcmp(s.name,{'primary','complex'}))
  ref=jsondecode(fileread(fullfile(repo,'soft_refit_20261023','results',s.name,'reference.json')));
  idx=strcmp({ref.datasets.name},s.condition);T=ref.datasets(idx).target;
  if strcmp(s.name,'primary'),weak=jsondecode(fileread(fullfile(repo,'soft_refit_20261023','inputs','primary.json')));end
 else
  tr=load(fullfile(repo,'bae_q_20261024','truth','independent',[s.name '_truth.mat']));T=tr.tr(strcmp({tr.tr.name},s.condition)).target;
 end
 live=any(s.D~=0,1);e=s.FullCalRefit-T;es=s.QScanRefit-T;
 relative_all=norm(e,'fro')/norm(T,'fro');scan_all=norm(es,'fro')/norm(T,'fro');
 relative_old=norm(e(:,live),'fro')/norm(T(:,live),'fro');scan_old=norm(es(:,live),'fro')/norm(T(:,live),'fro');
 weak_base=NaN;weak_scan=NaN;
 if ~isempty(weak)
  weak_base=score_weak(s.FullCalRefit,weak.reflectivity,weak.event_indices_zero_based+1,[2 4 6 8],s.cm.fpeak,s.cm.dt);
  weak_scan=score_weak(s.QScanRefit,weak.reflectivity,weak.event_indices_zero_based+1,[2 4 6 8],s.cm.fpeak,s.cm.dt);
 end
 ix=find(P.model==s.name&P.condition==s.condition);
 for ii=reshape(ix,1,[])
  m=s.cm.t>=P.start_s(ii)-1e-9&s.cm.t<P.end_s(ii)-1e-9;j=P.trace(ii);
  P.base_mse(ii)=mean(e(m,j).^2);P.scan_mse(ii)=mean(es(m,j).^2);P.error_rms(ii)=sqrt(P.base_mse(ii));
  P.target_rms(ii)=sqrt(mean(T(m,j).^2));
 end
 summary(end+1,:)={s.name,s.condition,relative_all,scan_all,relative_old,scan_old,scan_old/relative_old, ...
  weak_base,weak_scan,s.scaninfo.selected_factor}; %#ok<AGROW>
end
writetable(P,fullfile(root,'analysis','scored_windows.csv'));
writetable(cell2table(summary,'VariableNames',{'model','condition','base_all','scan_all','base_old_live', ...
 'scan_old_live','q0_error_ratio','weak_base','weak_scan','selected_Q_factor'}),fullfile(root,'analysis','case_scores.csv'));
fprintf('DONE scoring %d windows\n',height(P));
end
