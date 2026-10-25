function P=cb_setup()
% Paths, constants and the ten scenes of PLAN_ZH.md. Every input is read from cache/inputs,
% a byte-verified copy of Git objects (scripts/snapshot_inputs.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
P.bae=fullfile(P.in,'bae_q_20261024');P.r3=fullfile(P.in,'round3_20261024');P.bv=fullfile(P.in,'qscan_bias_variance_20261024');
addpath(fullfile(P.bae,'core'),fullfile(P.bae,'source'),fullfile(P.r3,'core'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid, verbatim
P.threads=4;   % round 3 ran with 4; the coherence step is not bit-stable across thread counts (scripts/diag_threads.m)
P.coherence=struct('half_width',2,'dip_periods',.25);
P.rules=struct('Standard',struct('lambda_min',.5,'lambda_max',16));   % round-3 "Standard"
P.methods={'FullCalRefit','Standard','QScanRefit','QScan_Standard','SoftGroupRefit','QScan_Soft'};   % A B C D + two extras
% Scenes: model, SNR, Q scale, replicates, seed level (seed = 1e9 + 1e6*model + 1e5*level + r).
m={'F1','F2','F1','F2','F1','F2','F1','F2','F1','F2'};snr={14,14,8,8,14,14,8,8,20,20};
u=1+(-30)/100;qs={1,1,1,1,u,u,u,u,1,1};   % Q error written as in run_combo.m: Q0*(1+e/100)
R={100,100,100,100,50,50,50,50,40,40};lev={1,1,2,2,1,1,2,2,3,3};
P.scenes=struct('id',num2cell(1:10),'model',m,'snr_db',snr,'q_scale',qs,'R',R,'seed_level',lev);
for s=1:numel(P.scenes)
 sc=P.scenes(s);P.scenes(s).name=sprintf('S%02d_%s_Q%03d_%02ddB',s,sc.model,round(100*sc.q_scale),sc.snr_db);
 P.scenes(s).seed_base=1e9+1e6*str2double(sc.model(2))+1e5*sc.seed_level;
end
end
