% Development (seeds 901-902): selections for the known wavelet and W3 with anchors v_a = 0 and 0.01.
addpath('orig_wang2021'); L=orig_lib();
for Qt=[50 100]
  for s=[901 902]
    d=gen_data(L,Qt,s); tic; S=scan_q(d,L,L.f0);
    fprintf('Qt %d s %d known: Q1 %.1f Q2 %.1f (%.0fs)\n',Qt,s,1/S.v1,1/S.v2,toc);
    fprintf('   E1 %s\n',mat2str(round(S.E1,2)));
    win=L.windows.W3; ts=mean(win); fp=peak_freq(d,L,win);
    for va=[0 0.01]
      fw=fw_anchor(fp,ts,va); S=scan_q(d,L,fw);
      fprintf('  W3 fp %d va %.2f fw %.1f: Q1 %.1f Q2 %.1f  v1 %.4f v2 %.4f\n',fp,va,fw,1/S.v1,1/S.v2,S.v1,S.v2);
    end
  end
end
