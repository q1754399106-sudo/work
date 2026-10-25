function run_reps(s,first,last)
% Stage 1 worker: replicates first..last of scene s. Each replicate draws a new noise field
% from the archived generator, re-estimates sigma, and runs all estimators on that same gather
% (paired design). Outputs go to cache/outputs/<scene>/ as single; a replicate whose file exists
% is skipped, so an interrupted worker can be relaunched. Truth is used only for the log line.
P=cb_setup();maxNumCompThreads(P.threads);sc=load_scene(P,s);assert(first>=1&&last<=sc.R);
dir_=fullfile(P.root,'cache','outputs',sc.name);if ~isfolder(dir_),mkdir(dir_);end
ET=sum(sc.T(:).^2);started=tic;
for r=first:last
 out=fullfile(dir_,sprintf('r%03d.mat',r));if isfile(out),continue;end
 seed=sc.seed_base+r;D=sc.D0+noise_replica(sc.D0,sc.snr_db,seed);sig=estimate_sigmas(D,sc.cm.t);
 [G,info]=combo_outputs(P,D,sig,sc.cm,sc.Qused,sc.offset);
 for m=P.methods,G.(m{1})=single(G.(m{1}));end
 save([out '.part'],'G','info','seed','sig','-v7');movefile([out '.part'],out);
 rel=@(m)sqrt(sum((double(G.(m))-sc.T).^2,'all')/ET);
 fprintf('%s r=%d factor %.4f A %.4f B %.4f C %.4f D %.4f flags %d (%.0fs, total %.0fs)\n',sc.name,r,info.qscan_factor, ...
  rel('FullCalRefit'),rel('Standard'),rel('QScanRefit'),rel('QScan_Standard'),info.optimization_flags,info.seconds,toc(started));
end
fprintf('DONE %s %d-%d %.0fs\n',sc.name,first,last,toc(started));
end
