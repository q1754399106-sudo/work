% Test data (seeds 1201-1206), once per (Q, seed): first 24 traces and the W3 dominant frequency from all 384. Deterministic; not committed.
addpath('../../qscan_orig_solver_20261106/scripts/orig_wang2021'); addpath('../../qscan_orig_solver_20261106/scripts');
L=orig_lib(); w=L.windows.W3;
for s=1201:1206
  for Qt=[50 70 100 140]
    d=gen_data(L,Qt,s); fp=peak_freq(d,L,w); dd=d(:,1:L.nscan);
    save(sprintf('../cache/d_Q%d_s%d.mat',Qt,s),'dd','fp','-v7'); fprintf('Q%d s%d fp %g\n',Qt,s,fp);
  end
end
