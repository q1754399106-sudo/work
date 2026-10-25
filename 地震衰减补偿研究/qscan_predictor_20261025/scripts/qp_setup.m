function P=qp_setup()
% Paths, constants, the ten development models and the three scenes of PLAN_ZH.md (Octave).
% Inputs come from cache/inputs, a byte-verified copy of Git objects (snapshot_inputs.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
addpath(fullfile(fileparts(P.root),'octave_compat_20261025'),'-begin');warning('off','all');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'));
P.grid=exp(linspace(log(.6),log(2.2),13));P.R=6;
names={'G1','G2','G3','G4','G5','G6','G7','G8','F1','F2'};
P.models=struct('name',names,'index',num2cell([1:8 1 2]),'family',[repmat({'G'},1,8) {'F','F'}]);
% Scenes: Q0 14 dB, Q0 8 dB, Q0 x 0.7 14 dB (same noise as Q0 14 dB in both families).
P.scenes=struct('tag',{'Q100_14dB','Q100_08dB','Q070_14dB'},'snr_db',{14,8,14},'q_scale',{1,1,.7},'level',{1,2,1});
P.seed_G=@(m,level,r) 1300000000+10000000*m+100000*level+r;   % qscan_offset_20261025/scripts/os_setup.m
P.seed_F=@(m,level,r) 1e9+1e6*m+1e5*level+r;                 % combo_bias_variance_20261024/scripts/cb_setup.m
end
