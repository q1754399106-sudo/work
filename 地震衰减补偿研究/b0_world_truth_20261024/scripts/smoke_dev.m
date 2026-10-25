function smoke_dev()
% Development smoke test: diagnose2 on one F2/Moderate gather (id 15) with timing.
P=bw_setup();maxNumCompThreads(P.threads);
M=readtable(fullfile(P.in,'qscan_bias_diagnostic_20261024','predictions','manifest.csv'),'TextType','string');
g=M(M.id==15,:);[D,cm,Qused,offset]=observed_gather(P,g);r=diagnose2(P,D,cm,Qused,offset,g.boot_seed);disp(r);
end
