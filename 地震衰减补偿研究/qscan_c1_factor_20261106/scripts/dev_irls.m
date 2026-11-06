% Development only (seeds 901-902, as the earlier experiment's development seeds): IRLS with eq. 17, chengD, C2; exports data and A for the Python port.
addpath('../../qscan_orig_solver_20261106/scripts/orig_wang2021'); addpath('../../qscan_orig_solver_20261106/scripts');
L=orig_lib(); mkdir('../cache');
[~,~,A]=qcdn_core(zeros(L.n1,1),L.dt,24,1,1,70); save('../cache/dev_A_fw24_Q70.mat','A','-v7');
for s=[901 902]
 for Qt=[70 100 140]
  d=gen_data(L,Qt,s); w=L.windows.W3; fp=peak_freq(d,L,w); dd=d(:,1:L.nscan); save(sprintf('../cache/dev_d_Q%d_s%d.mat',Qt,s),'dd','fp','-v7');
  conds={'known',30;'W3_0',fw_anchor(fp,mean(w),0);'W3_1',fw_anchor(fp,mean(w),0.01);'F0.8',24;'F1.25',37.5};
  for c=1:size(conds,1)
   t0=tic; S=scan_q2(d,L,conds{c,2});
   fprintf('Q%d s%d %s fw %.1f | e17 Q %.0f cross %d ratio(Qsel/(Q fw/f0)) %.2f | chengD Q %.0f | C2 Q %.0f | %.0fs\n',Qt,s,conds{c,1},conds{c,2},1/S.v1,S.cross,(1/S.v1)/(Qt*conds{c,2}/30),1/S.vd,1/S.v2,toc(t0));
  end
 end
end
