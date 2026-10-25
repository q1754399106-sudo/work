function run_truth_replicates(set,model,condition)
% Scoring side, run only after the predictions are sealed: expected Delta MSE of QScan versus
% FullCalRefit at the public Q for one confirmation scene, from 12 noise replicates (the two
% diagnosed replicate gathers are replicates 1-2 of this set).
P=bd_setup();maxNumCompThreads(P.threads);assert(isfile(fullfile(P.root,'predictions','SEALED_SHA256.json')));
o=load_obs(P,set,model,condition);[D0,T]=regen_clean(P,set,model);live=any(o.D~=0,1);
snr=14;k=1;if strcmp(condition,'Strong'),snr=8;k=2;end
s=find(strcmp(set,{'bae','qc','r3'}));m=str2double(model(2));rows={};
for r=1:12
 seed=1100000000+10000000*s+1000000*m+100000*k+r;D=D0+noise_replica(D0,snr,seed);sig=estimate_sigmas(D,o.cm.t);
 sel=qscan_select(D,sig,o.cm,o.Q0,o.offset);c=P.grid(sel.best);
 x1=full_refit(D,sig,o.cm,o.Q0,o.offset);xs=full_refit(D,sig,o.cm,o.Q0*c,o.offset);ET=sum(T(:,live).^2,'all');
 rows(end+1,:)={set,model,condition,r,seed,c,sum((x1(:,live)-T(:,live)).^2,'all')/ET,sum((xs(:,live)-T(:,live)).^2,'all')/ET}; %#ok<AGROW>
end
writetable(cell2table(rows,'VariableNames',{'set','model','condition','replicate','seed','qscan_factor','mse_base','mse_qscan'}), ...
 fullfile(P.root,'analysis',sprintf('truth_replicates_%s_%s_%s.csv',set,model,condition)));
fprintf('DONE truth replicates %s %s %s\n',set,model,condition);
end
