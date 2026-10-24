function posthoc_b0()
% POST HOC (after E1 failed; not a registered test). On the development Q0 replicate gathers:
% (1) bootstrap noise model versus the true noise D - D0; (2) B0 rerun with the true per-trace
% noise level and the true AR coefficient (oracle noise), same seeds and K; (3) signal energy
% of the FullCalRefit output relative to the target, to put B0 on the target's scale.
P=bd_setup();maxNumCompThreads(P.threads);M=readtable(fullfile(P.root,'predictions','manifest.csv'),'TextType','string');
Dg=readtable(fullfile(P.root,'predictions','diagnostics_w1.csv'));for w=2:4,Dg=[Dg;readtable(fullfile(P.root,'predictions',sprintf('diagnostics_w%d.csv',w)))];end %#ok<AGROW>
sel=M(M.set=="bae"&M.q_err==0&ismember(M.condition,["Moderate","Strong"]),:);rows={};
for i=1:height(sel)
 g=sel(i,:);[D,cm,Qused,offset]=observed_gather(P,g);[D0,T]=regen_clean(P,'bae',char(g.model));N=D-D0;
 [phi,sigb,tend]=noise_model(D,cm);m=cm.t>=.02&cm.t<=tend;sig_true=std(N,0,1);
 Nz=N-mean(N,1);phi_true=sum(Nz(2:end,:).*Nz(1:end-1,:),'all')/sum(Nz.^2,'all');
 L=load(fullfile(P.root,'cache','outputs',sprintf('g%03d.mat',g.id)));x1=double(L.x1);
 r=Dg(Dg.id==g.id,:);c=r.qscan_factor;seeds=g.boot_seed+(1:P.K);
 H=bootstrap_world(x1,Qused,Qused,Qused*c,cm,offset,phi_true,sig_true,seeds);
 rows(end+1,:)={g.id,char(g.model),char(g.condition),mean(sigb./sig_true),phi,phi_true,sum(x1(:).^2)/sum(T(:).^2), ...
  r.B0_rel_delta_mse,r.B0_delta_bias2,r.B0_delta_var,H.rel_delta_mse,H.bias2(2)-H.bias2(1),H.var(2)-H.var(1),nnz(m)}; %#ok<AGROW>
 fprintf('g%03d %s %s: sigma_hat/true %.3f phi %.3f (true %.3f) |x1|^2/|T|^2 %.3f | B0 %+.3f -> oracle-noise %+.3f (dB2 %+.4f dVar %+.4f)\n', ...
  rows{end,1},rows{end,2},rows{end,3},rows{end,4},rows{end,5},rows{end,6},rows{end,7},rows{end,8},rows{end,11},rows{end,12},rows{end,13});
end
writetable(cell2table(rows,'VariableNames',{'id','model','condition','sigma_hat_over_true','phi_hat','phi_true','output_energy_over_target', ...
 'B0_rel_delta_mse','B0_delta_bias2','B0_delta_var','oracle_noise_rel_delta_mse','oracle_noise_delta_bias2','oracle_noise_delta_var','window_samples'}), ...
 fullfile(P.root,'analysis','posthoc_b0_noise.csv'));
fprintf('DONE posthoc\n');
end
