function P=os_setup()
% Paths, constants, the eight development geologies G1-G8 and the three scenes of PLAN_ZH.md.
% Every input is read from cache/inputs, a byte-verified copy of Git objects (snapshot_inputs.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'),fullfile(P.in,'round3_20261024','core'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid, verbatim
P.factors=[1 P.grid];P.R=6;
P.threads=4;   % as combo_bias_variance_20261024 (the coherence step is not bit-stable across thread counts)
P.coherence=struct('half_width',2,'dip_periods',.25);
P.rules=struct('Standard',struct('lambda_min',.5,'lambda_max',16));   % round-3 "Standard"
% Development geologies: make_unseen.m / make_independent.m generator, new seeds and macro models
% (none equal to F1, F2 or U1-U4).
P.models=struct('name',{'G1','G2','G3','G4','G5','G6','G7','G8'}, ...
 'seed',{2026102601,2026102611,2026102621,2026102631,2026102641,2026102651,2026102661,2026102671}, ...
 'macro_vp',{[1850 2300 2750 3300],[1900 2350 2800 3350],[1780 2200 2820 3380],[1830 2380 2720 3260], ...
             [1870 2260 2790 3320],[1810 2340 2760 3290],[1760 2280 2840 3360],[1890 2310 2740 3240]}, ...
 'macro_h',{[240 300 350],[230 320 330],[250 290 360],[235 305 345],[255 315 325],[245 285 365],[225 335 335],[265 275 350]}, ...
 'macro_Q',{[50 70 95 130],[90 60 45 40],[40 80 50 110],[60 50 70 90],[30 45 60 80],[120 100 90 80],[70 40 100 60],[45 90 60 140]}, ...
 'carbonate',{false,true,false,true,false,true,false,true});
% Scenes: Q0 at 14 and 8 dB; Q0 x 0.7 at 14 dB on the same noise as the Q0 14 dB scene.
P.scenes=struct('tag',{'Q100_14dB','Q100_08dB','Q070_14dB'},'snr_db',{14,8,14},'q_scale',{1,1,.7},'level',{1,2,1});
P.seed=@(m,level,r) 1300000000+10000000*m+100000*level+r;
end
