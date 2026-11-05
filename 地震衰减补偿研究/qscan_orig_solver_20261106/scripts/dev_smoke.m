% Development (seed 901, not a test seed): are the criterion curves usable, how long does a scan take?
addpath('orig_wang2021'); L=orig_lib();
for Qt=[50 100]
    d=gen_data(L,Qt,901);
    for w=fieldnames(L.windows)'
        fprintf('Qt %d %s fp %d\n',Qt,w{1},peak_freq(d,L,L.windows.(w{1})));
    end
    tic; S=scan_q(d,L,L.f0); el=toc;
    fprintf('Qt %d known wavelet: Q(C1) %.1f Q(C2) %.1f  %.0fs\n',Qt,1/S.v1,1/S.v2,el);
    fprintf('  C1 %s\n  C2 %s\n',mat2str(round(S.C1,3)),mat2str(round(S.C2,1)));
end
