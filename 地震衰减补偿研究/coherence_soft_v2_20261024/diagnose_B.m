function diagnose_B()
% Stage B diagnosis (truth side): where does the remaining V3 error sit?
root=fileparts(mfilename('fullpath'));base=fileparts(root);v1=fullfile(base,'coherence_soft_20261023');
addpath(fullfile(root,'core'),fullfile(v1,'core'),fullfile(v1,'source'));
TT=load(fullfile(v1,'truth','V3_truth.mat'),'g');in=load(fullfile(v1,'inputs','validation','V3.mat'),'obs');
maxNumCompThreads(4);model=gen_complex_components(TT.g);
for k=1:2,dv=norm(model.datasets(k).D-in.obs(k).D,'fro')/norm(in.obs(k).D,'fro');fprintf('regenerated D relative difference %.3g\n',dv);assert(dv<1e-10,'Regenerated V3 differs.');end
best=struct('lambda_min',.5,'lambda_max',16,'tau',0,'kappa',0);rows={};
for e=[-30 0 30]
 for k=1:2
  nm=model.datasets(k).name;S=load(fullfile(root,'cache','scenes',sprintf('V3_Q%+d_%s.mat',e,nm)),'s');s=S.s;s.R0=full(s.R0);
  d=prune_dictionary(s.cm,450);s.W=d.W;A=apply_gain(s,s.c,best);G=load(fullfile(root,'cache','truth',sprintf('V3_Q%+d_%s.mat',e,nm)),'gathers');
  ds=model.datasets(k);sm=@(x)sqrt(movmean(x.^2,21,1));
  P=sm(ds.primary_att);Cc=sm(ds.coherent);Mm=sm(ds.multiple);Nn=sm(ds.noise);
  [~,cls]=max(cat(3,P,Cc,Mm,Nn),[],3);cls(:,~model.live)=0;names={'primary','coherent','multiple','noise'};
  for method={'FullCalRefit','SoftGroupRefit','best'}
   if strcmp(method{1},'best'),X=A;else,X=G.gathers.(method{1});end
   err=(X-ds.target).^2;tot=sum(err(:));
   for c=1:4,rows{end+1}=struct('Q',e,'condition',nm,'method',method{1},'class',names{c}, ...
     'sample_share',mean(cls(:)==c),'error_share',sum(err(cls==c))/tot,'error_energy',sum(err(cls==c))/sum(ds.target(:).^2)); %#ok<AGROW>
   end
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis','stage_B_diagnosis.csv'));disp(M);
end
