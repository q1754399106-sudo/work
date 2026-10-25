function P=cc_setup()
% Paths, constants and the 30 scenes of PLAN_ZH.md: F1, F2 (combo_bias_variance_20261024 stage 2)
% and G1-G8 (qscan_offset_20261025), each at Q0/14 dB, Q0/8 dB and Q0x0.7/14 dB, replicates 1-6.
% Model inputs come from cache/inputs (byte-verified Git objects, snapshot_inputs.py); the
% per-factor coherence Standard outputs come from the two archived sweep caches (hashes in
% cache/SNAPSHOT.json).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');P.repo=fileparts(P.root);
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid, verbatim (the sweeps hold outputs at 1 and these 13)
P.factors=[1 P.grid];P.R=6;P.threads=4;
P.criteria={'A0','A1','A2','A3','C1','C2','C3','C4'};
tags={'Q100_14dB','Q100_08dB','Q070_14dB'};snr=[14 8 14];qs=[1 1 .7];
P.scenes=struct('name',{},'model',{},'family',{},'tag',{},'snr_db',{},'q_scale',{},'primary',{},'sweep',{});
combo={'S01_F1_Q100_14dB','S03_F1_Q100_08dB','S05_F1_Q070_14dB';'S02_F2_Q100_14dB','S04_F2_Q100_08dB','S06_F2_Q070_14dB'};
for m=1:2
 for s=1:3
  P.scenes(end+1)=struct('name',combo{m,s},'model',sprintf('F%d',m),'family','F','tag',tags{s},'snr_db',snr(s),'q_scale',qs(s), ...
   'primary',snr(s)==14,'sweep',fullfile(P.repo,'combo_bias_variance_20261024','cache','sweep',combo{m,s}));
 end
end
for m=1:8
 for s=1:3
  name=sprintf('G%d_%s',m,tags{s});
  P.scenes(end+1)=struct('name',name,'model',sprintf('G%d',m),'family','G','tag',tags{s},'snr_db',snr(s),'q_scale',qs(s), ...
   'primary',snr(s)==14,'sweep',fullfile(P.repo,'qscan_offset_20261025','cache','sweep',name));
 end
end
end
