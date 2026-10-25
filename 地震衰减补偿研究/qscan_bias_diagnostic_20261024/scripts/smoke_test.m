function smoke_test()
% Development smoke test on the first manifest gather (F1/Moderate replicate 1): runs diagnose
% once, prints its statistics and timing; writes nothing into predictions/.
P=bd_setup();maxNumCompThreads(P.threads);M=readtable(fullfile(P.root,'predictions','manifest.csv'),'TextType','string');
g=M(1,:);[D,cm,Qused,offset]=observed_gather(P,g);tic;r=diagnose(P,D,cm,Qused,offset,g.boot_seed);r.seconds=toc;disp(r);
end
