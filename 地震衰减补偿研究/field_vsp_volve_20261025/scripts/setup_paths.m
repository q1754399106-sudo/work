function P=setup_paths()
% Paths and constants of PLAN_ZH.md. Frozen code is read from cache/inputs (scripts/snapshot_inputs.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
addpath(fullfile(P.root,'scripts','core'),fullfile(P.in,'bae_q_20261024','core'),fullfile(P.in,'bae_q_20261024','source'), ...
 fullfile(P.in,'round3_20261024','core'),fullfile(P.in,'combo_bias_variance_20261024','scripts'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid of combo_outputs/qscan_select, verbatim
P.threads=4;
P.coherence=struct('half_width',2,'dip_periods',.25);
P.rules=struct('Standard',struct('lambda_min',.5,'lambda_max',16));
P.methods={'FullCalRefit','Standard','QScanRefit','QScan_Standard','SoftGroupRefit','QScan_Soft'};
P.runs=struct('id',{'V0','Vm50','Vm30','Vp30','Vp50','R0'},'model',{'VSP','VSP','VSP','VSP','VSP','READ'}, ...
 'e',{0,-50,-30,30,50,0});
P.levels=191;P.seed=20261025;
end
