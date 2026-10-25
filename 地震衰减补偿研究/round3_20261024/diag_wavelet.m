function diag_wavelet()
% Interim diagnosis 1: how much of the complex-model error is the zero-phase
% wavelet assumption? Oracle (true mean phase 24 deg + ghost) vs Ricker, plus a
% data-only phase scan with the QScan-style lasso objective (no truth used).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
in=load(fullfile(root,'dev','V5.mat'));T=load(fullfile(root,'dev','V5_truth.mat'),'tr');cm=in.config;cm.t=cm.t(:);
d0=prune_dictionary(cm,450);dz=wavelet_dictionary(cm,450,0);
fprintf('zero-phase wavelet_dictionary vs prune_dictionary: %.3g\n',norm(dz.B-d0.B,'fro')/norm(d0.B,'fro'));
ghost=[-.16 .012;.04 .025];
variants={'Ricker',d0;'Phase24',wavelet_dictionary(cm,450,24);'Phase24Ghost',wavelet_dictionary(cm,450,24,ghost)};
rows={};
for k=1:2
 D=in.obs(k).D;sig=reshape(in.obs(k).sigma,1,[]);tg=T.tr(k).target;
 for v=1:size(variants,1)
  d=variants{v,2};S=zeros(size(D));
  for b=1:9,ids=(b-1)*11+(1:11);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,11),sig(ids));S(:,ids)=o.gather;end
  e=norm(S-tg,'fro')/norm(tg,'fro');rows{end+1}=sprintf('%s %s FullCalRefit error %.4f',in.obs(k).name,variants{v,1},e); %#ok<AGROW>
  fprintf('%s\n',rows{end});
 end
 % Data-only phase scan: gather-summed lasso objective, as in qscan_refit.
 n=size(D,1);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sig;phases=0:8:56;J=zeros(size(phases));
 for p=1:numel(phases)
  d=wavelet_dictionary(cm,450,phases(p));X=refit_lasso(d.Gram,d.B'*D(:,live),lambda(live),d.L);res=D(:,live)-d.B*X;
  J(p)=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sig(live).^2);
 end
 [~,b]=min(J);fprintf('%s phase scan J(%s) = %s -> %d deg\n',in.obs(k).name,mat2str(phases),mat2str(round(J)),phases(b));
end
end
