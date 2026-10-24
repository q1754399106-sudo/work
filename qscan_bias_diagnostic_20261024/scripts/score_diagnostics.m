function score_diagnostics()
% Scoring side, after sealing: realized errors of both estimators on every diagnosed gather.
P=bd_setup();assert(isfile(fullfile(P.root,'predictions','SEALED_SHA256.json')));
M=readtable(fullfile(P.root,'predictions','manifest.csv'),'TextType','string');rows={};targets=containers.Map();
for i=1:height(M)
 g=M(i,:);o=load_obs(P,char(g.set),char(g.model),'Moderate');key=char(g.set+"/"+g.model);
 if ~isKey(targets,key),[~,targets(key)]=regen_clean(P,char(g.set),char(g.model));end
 T=targets(key);
 live=any(o.D~=0,1);L=load(fullfile(P.root,'cache','outputs',sprintf('g%03d.mat',g.id)));ET=sum(T(:,live).^2,'all');
 rows(end+1,:)={g.id,sum((double(L.x1(:,live))-T(:,live)).^2,'all')/ET,sum((double(L.xs(:,live))-T(:,live)).^2,'all')/ET}; %#ok<AGROW>
end
writetable(cell2table(rows,'VariableNames',{'id','mse_base','mse_qscan'}),fullfile(P.root,'analysis','realized_scores.csv'));
fprintf('DONE scoring %d gathers\n',height(M));
end
