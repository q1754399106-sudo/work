function run_clean()
% Stage 1: deterministic counterpart of every estimator. Noise-free gather D0 with the
% archived sigma, i.e. the same lambda as the archived run but no noise.
P=bv_setup();maxNumCompThreads(P.threads);rows={};sel_rows={};nf=numel(P.factors);
for s=1:numel(P.scenes)
 sc=load_scene(P,s);live=any(sc.D~=0,1);X=zeros([size(sc.D0) nf]);
 sel=qscan_select(sc.D0,sc.sigma,sc.cm,sc.Q0,sc.offset);
 for f=1:nf
  [X(:,:,f),info]=full_refit(sc.D0,sc.sigma,sc.cm,sc.Q0*P.factors(f),sc.offset);
  rel=norm(X(:,live,f)-sc.T(:,live),'fro')/norm(sc.T(:,live),'fro');
  rows(end+1,:)={sc.model,sc.condition,f,P.factors(f),rel,mean(info.support_size),info.optimization_flags,info.max_normal_residual}; %#ok<AGROW>
 end
 for q=1:numel(P.grid)
  sel_rows(end+1,:)={sc.model,sc.condition,q,P.grid(q),sel.objective(q),sel.misfit(q),sel.penalty(q),sel.lasso_active(q),q==sel.best}; %#ok<AGROW>
 end
 X=single(X);save(fullfile(P.root,'results',sprintf('clean_%s_%s.mat',sc.model,sc.condition)),'X','sel','-v7');
 fprintf('%s %s clean: selected %.4f, rel(c=1) %.4f, rel(selected) %.4f\n',sc.model,sc.condition,sel.selected_factor, ...
  rows{end-nf+1,5},rows{end-nf+1+sel.best,5});
end
writetable(cell2table(rows,'VariableNames',{'model','condition','slot','factor','relative_error','mean_support', ...
 'optimization_flags','max_normal_residual'}),fullfile(P.root,'analysis','clean_sweep.csv'));
writetable(cell2table(sel_rows,'VariableNames',{'model','condition','grid_index','factor','objective','misfit','penalty', ...
 'lasso_active_per_trace','selected'}),fullfile(P.root,'analysis','clean_qscan_objective.csv'));
fprintf('DONE stage 1\n');
end
