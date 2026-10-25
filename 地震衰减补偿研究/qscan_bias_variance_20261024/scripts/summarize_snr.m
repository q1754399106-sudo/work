function summarize_snr()
% Addendum A: bias^2/variance split for the high-SNR probe (same corrections as summarize.m).
P=bv_setup();rows={};
for model_number=1:2
 s=find(strcmp({P.scenes.model},sprintf('F%d',model_number))&strcmp({P.scenes.condition},'Moderate'));
 sc=load_scene(P,s);live=any(sc.D~=0,1);T=sc.T;ET=sum(T(:,live).^2,'all');
 for snr=[20 26]
  files=dir(fullfile(P.root,'results',sprintf('snr_%s_%02ddB_*.mat',sc.model,snr)));
  M1=0;M2=0;REL=[];rep=[];best=[];
  for f=reshape(files,1,[])
   L=load(fullfile(f.folder,f.name));M1=M1+L.M1;M2=M2+L.M2;REL=[REL;L.REL];rep=[rep;L.replicate];best=[best;L.best]; %#ok<AGROW>
   labels=L.labels;
  end
  assert(isequal(sort(rep),(1:40)'),'Probe replicates missing or duplicated.');R=numel(rep);
  Xbar=M1/R;V=max(M2/R-Xbar.^2,0)*R/(R-1);
  for k=1:numel(labels)
   va=sum(V(:,live,k),'all');b2=sum((Xbar(:,live,k)-T(:,live)).^2,'all')-va/R;mse=b2+va;
   assert(abs(mean(REL(:,k).^2)*ET-mse)/mse<1e-9,'MSE identity failed.');
   ratio=REL(:,k)./REL(:,1);
   rows(end+1,:)={sc.model,snr,labels{k},R,mse/ET,b2/ET,va/ET,mean(REL(:,k)),mean(ratio>1),median(ratio), ...
    mean(best==8),numel(unique(best))}; %#ok<AGROW>
  end
  fprintf('%s %d dB: base %.4f fixed1.28 %.4f QScan %.4f (mse rel)\n',sc.model,snr,rows{end-2,5},rows{end-1,5},rows{end,5});
 end
end
writetable(cell2table(rows,'VariableNames',{'model','snr_db','estimator','R','mse_rel','bias2_rel','var_rel','mean_rel_error', ...
 'fraction_worse_than_base','median_ratio_to_base','fraction_selected_1p2803','distinct_selected'}),fullfile(P.root,'analysis','snr_probe.csv'));
fprintf('DONE snr summary\n');
end
