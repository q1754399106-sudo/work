function run_field(line)
% All registered methods on the prepared window, sediment-Q errors e (PLAN_ZH.md).
% QScan additionally at e=+40,+60 (starts 1.4 and 1.6 x Q_SR, registered E4-2).
root=fileparts(mfilename('fullpath'));bae=fullfile(fileparts(root),'bae_q_20261024');
addpath(fullfile(root,'core'),fullfile(bae,'core'),fullfile(bae,'source'));maxNumCompThreads(4);
in=load(fullfile(root,'inputs',[line '.mat']));D=in.D;cm=in.cm;sigma=in.sigma;refQ=cm.Q;
opts=struct('s_grid',.10,'M',32,'rho',.7,'seed',20261024,'ensemble_nodes',7,'lean',true);
outdir=fullfile(root,'results');if ~isfolder(outdir),mkdir(outdir);end;started=tic;
for e=[-50 -30 0 30 50 40 60]
 file=fullfile(outdir,sprintf('%s_Q%+d.mat',line,e));if isfile(file),continue;end
 used=cm;used.Q=refQ*(1+e/100);g=struct('Input',D);
 [S,qinfo]=qscan_refit(D,sigma,used,0);g.QScanRefit=S;info=struct('qscan',qinfo,'used_Q',used.Q);
 if ismember(e,[-50 -30 0 30 50])
  [gs,diag]=run_scene(cm,D,sigma,refQ,e,opts);g.FullCalRefit=gs.FullCalRefit;g.SoftGroupRefit=gs.SoftGroupRefit;
  g.QEnsembleMean_s010=gs.QEnsembleMean_s010;info.scene=diag;
  d=prune_dictionary(used,0);sbl=zeros(size(D));
  for first=1:11:size(D,2),ids=first:min(first+10,size(D,2));o=sbl_tvd(D(:,ids),d,sigma(ids));sbl(:,ids)=o.gather;end
  g.SBL_TVD=sbl;g.WangIQF_G10=wang_inverse_q(D,used,0,10);
 end
 save(file,'g','info','-v7');
 fprintf('%s e=%+d QScan factor %.3f implied Q %.1f elapsed %.1fs\n',line,e,qinfo.selected_factor,qinfo.applied_Q(2),toc(started));
end
fprintf('DONE %s\n',line);
end
