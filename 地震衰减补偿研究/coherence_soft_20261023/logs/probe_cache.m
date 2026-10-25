% Development probe: structure of cached soft-refit checkpoints (read-only).
root='D:/ChatGPT/研究/soft_refit_20261023';
c=load(fullfile(root,'cache','primary_Q+0_Strong_05.mat'));
disp(fieldnames(c.initial)');disp(fieldnames(c.out)');
in=jsondecode(fileread(fullfile(root,'inputs','primary.json')));
ev=in.event_indices_zero_based+1;
ng=cellfun(@numel,c.out.group_gain);fprintf('groups per trace: %s\n',mat2str(ng));
cm=in.config;gap=max(1,round(.5/(cm.fpeak*cm.dt)));fprintf('gap=%d\n',gap);
for q=1:11
 j=44+q;s=c.initial.support{q}(:);g=cumsum([1;diff(s)>gap]);
 cen=accumarray(g,s,[],@mean);gain=c.out.group_gain{q};
 isev=arrayfun(@(x)any(abs(ev(j,:)-x)<=gap),cen);
 fprintf('trace %d: %d groups, %d near events; gain(event)=%s gain(other)=%s\n',j,numel(cen),nnz(isev), ...
   mat2str(round(gain(isev)',2)),mat2str(round(gain(~isev)',2)));
end
fprintf('weak event amplitudes at trace 50: %s\n',mat2str(in.reflectivity(sub2ind(size(in.reflectivity),ev(50,:),repmat(50,1,9))),3));
