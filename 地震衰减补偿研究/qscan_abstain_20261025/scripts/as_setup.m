function P=as_setup()
% Paths, constants, the eight new confirmation geologies K1-K8 and the scenes of PLAN_ZH.md (Octave).
% Inputs come from cache/inputs, a byte-verified copy of Git objects (snapshot_inputs.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
addpath(fullfile(fileparts(P.root),'octave_compat_v2_20261025'),'-begin');warning('off','all');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'),fullfile(P.in,'round3_20261024','core'));
P.grid=exp(linspace(log(.6),log(2.2),13));P.R=6;
P.coherence=struct('half_width',2,'dip_periods',.25);
P.rules=struct('Standard',struct('lambda_min',.5,'lambda_max',16));   % round-3 "Standard"
P.rule=struct('feature','Z1','threshold',72.65656565656568,'shift',2);   % frozen in PLAN_ZH.md
% Confirmation geologies: make_dev.m generator, new seeds (none used by F, G, H, U) and macro models.
P.models=struct('name',{'K1','K2','K3','K4','K5','K6','K7','K8'}, ...
 'seed',{2026102801,2026102811,2026102821,2026102831,2026102841,2026102851,2026102861,2026102871}, ...
 'macro_vp',{[1840 2290 2770 3310],[1880 2330 2790 3340],[1790 2230 2830 3370],[1820 2370 2740 3280], ...
             [1860 2270 2800 3300],[1800 2350 2770 3310],[1770 2290 2850 3350],[1900 2320 2730 3250]}, ...
 'macro_h',{[235 310 340],[245 300 355],[240 325 330],[255 290 350],[230 320 345],[250 295 355],[245 330 325],[260 280 360]}, ...
 'macro_Q',{[55 75 100 125],[85 65 50 45],[45 85 55 105],[65 55 75 95],[35 50 65 85],[110 95 85 75],[75 45 95 65],[50 85 65 130]}, ...
 'carbonate',{false,true,false,true,false,true,false,true});
% Scenes: Q0 14 dB and Q0 x 0.7 14 dB (primary, same noise); Q0 8 dB and Q0 20 dB (secondary).
P.scenes=struct('tag',{'Q100_14dB','Q070_14dB','Q100_08dB','Q100_20dB'},'snr_db',{14,14,8,20},'q_scale',{1,.7,1,1},'level',{1,1,2,3});
P.seed=@(m,level,r) 1500000000+10000000*m+100000*level+r;
end
