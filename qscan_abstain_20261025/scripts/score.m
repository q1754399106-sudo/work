function score()
% Scoring (after all decisions are sealed): per gather, MSE/||x*||^2 of B, D, S2 and of the rule.
P=as_setup();assert(isfile(fullfile(P.root,'decisions','SEALED_SHA256.json')));
fid=fopen(fullfile(P.root,'analysis','gather_scores.csv'),'w');fprintf(fid,'model,scene,replicate,b,g,Z1,mse_B,mse_D,mse_S2,mse_rule\n');
for m=1:numel(P.models)
 T=load(fullfile(P.root,'truth',[P.models(m).name '_truth.mat']),'target');T=T.target;ET=sum(T(:).^2);
 for s=1:numel(P.scenes)
  name=[P.models(m).name '_' P.scenes(s).tag];
  for r=1:P.R
   L=load(fullfile(P.root,'cache','out',name,sprintf('r%03d.mat',r)));e=zeros(1,3);
   for k=1:3,X=double(L.STD(:,:,k));e(k)=sum((X(:)-T(:)).^2)/ET;end
   if L.g==L.b,er=e(2);else,er=e(3);end
   fprintf(fid,'%s,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n',P.models(m).name,P.scenes(s).tag,r,L.b,L.g,L.Z1,e(1),e(2),e(3),er);
  end
 end
end
fclose(fid);disp('DONE score');
end
