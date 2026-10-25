function run_diagnostics(worker,workers)
% Stage 2 worker: manifest rows with mod(id-1,workers)==worker-1. Writes truth-free predictions
% (one CSV per worker) and the two outputs per gather for later scoring. No truth is read.
P=bd_setup();maxNumCompThreads(P.threads);M=readtable(fullfile(P.root,'predictions','manifest.csv'),'TextType','string');
M=M(mod(M.id-1,workers)==worker-1,:);out=fullfile(P.root,'predictions');outs=fullfile(P.root,'cache','outputs');
if ~isfolder(outs),mkdir(outs);end
csv=fullfile(out,sprintf('diagnostics_w%d.csv',worker));assert(~isfile(csv),'Worker output exists; refusing to overwrite.');
rows=[];started=tic;
for i=1:height(M)
 g=M(i,:);t0=tic;[D,cm,Qused,offset]=observed_gather(P,g);
 [r,x1,xs]=diagnose(P,D,cm,Qused,offset,g.boot_seed);x1=single(x1);xs=single(xs); %#ok<NASGU>
 save(fullfile(outs,sprintf('g%03d.mat',g.id)),'x1','xs','-v7');
 r.id=g.id;r.seconds=toc(t0);rows=[rows;struct2table(r)]; %#ok<AGROW>
 writetable(rows,csv);
 fprintf('g%03d %s %s %s q%+d %s: c=%.4f B0 %+.3f gain %.4f cv %+.4f rs %+.4f (%.0fs, total %.0fs)\n',g.id,g.set,g.model, ...
  g.condition,g.q_err,g.source,r.qscan_factor,r.B0_rel_delta_mse,r.gate_gain,r.cv_rel_delta,r.rs_delta,r.seconds,toc(started));
end
fprintf('DONE worker %d: %d gathers %.0fs\n',worker,height(M),toc(started));
end
