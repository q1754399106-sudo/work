% Demo: both modes on one bundled observation (independent Kennett forward F4,
% Strong noise) with the calibrated Q deliberately 30% too low.
here=fileparts(mfilename('fullpath'));root=fileparts(here);addpath(root);
in=load(fullfile(root,'tests','data','F4.mat'));ob=in.obs(strcmp({in.obs.name},'Strong'));
cm=in.config;cm.Q=cm.Q*0.7;
a=q_compensate(ob.D,cm,'sigma',ob.sigma);                    % Q as given
b=q_compensate(ob.D,cm,'sigma',ob.sigma,'mode','qscan');     % data-driven common Q factor
fprintf('QScan factor %.3f -> applied Q %s (given %s)\n',b.qscan_factor,mat2str(b.applied_Q,4),mat2str(cm.Q,4));
fprintf('fallback groups: %d of %d\n',sum(cellfun(@nnz,b.fallback)),sum(cellfun(@numel,b.fallback)));
if usejava('desktop')||usejava('jvm')
 t=cm.t;c=prctile(abs(b.gather(:)),99.5);figure('Position',[50 50 1300 450]);
 subplot(1,3,1);imagesc(1:size(ob.D,2),t,ob.D,[-1 1]*prctile(abs(ob.D(:)),99.5));title('observed');ylabel('time (s)');xlabel('trace');
 subplot(1,3,2);imagesc(1:size(ob.D,2),t,a.gather,[-c c]);title('standard (Q 30% low)');xlabel('trace');
 subplot(1,3,3);imagesc(1:size(ob.D,2),t,b.gather,[-c c]);title(sprintf('qscan + standard (c = %.2f)',b.qscan_factor));xlabel('trace');colormap(gray);
end
