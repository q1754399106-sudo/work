function P=fd_setup()
% fd2d_confirm3_20261114: fd2d_confirm_20261109/scripts/fd_setup.m with new models H13-H18 only (UQW-S confirmation).
% Paths, constants, the six 2-D development geologies H1-H6 and the three scenes of PLAN_ZH.md.
% Every upstream input is read from cache/inputs, a byte-verified copy of Git objects (snapshot_inputs.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'),fullfile(P.in,'round3_20261024','core'));
% Sweep constants: verbatim from qscan_offset_20261025/scripts/os_setup.m.
P.grid=exp(linspace(log(.6),log(2.2),13));P.factors=[1 P.grid];P.R=6;P.threads=4;
P.coherence=struct('half_width',2,'dip_periods',.25);
P.rules=struct('Standard',struct('lambda_min',.5,'lambda_max',16));
% Solver-side constants of every Kennett data set (F1-F4, G1-G8).
P.cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32,'lambda_x',.001,'lambda_fixed',.001);P.cm.t=(0:P.cm.nt-1)'*P.cm.dt;
% 2-D finite-difference forward (fd2d_visco): grid, time step, window around each zero-offset shot.
P.fd=struct('dx',4,'dt',4e-4,'npml',25,'half_window',500,'zmax',2950,'t0',.06,'nt',4700, ...
 'trace_dx',12,'ntraces',99,'r0',500,'sub',.25);
% Development geologies: fine-layer generator of make_independent.m / make_dev.m (new seeds and
% macro models, none equal to F1-F4, G1-G8 or U1-U4), laid out in 2-D with a vertical fault.
P.models=struct('name',{'H13','H14','H15','H16','H17','H18'}, ...
 'seed',{2026111401,2026111411,2026111421,2026111431,2026111441,2026111451}, ...
 'macro_vp',{[1860 2290 2760 3310],[1850 2330 2790 3340],[1820 2270 2810 3350],[1870 2350 2750 3300],[1830 2300 2780 3320],[1840 2320 2770 3280]}, ...
 'macro_h',{[245 295 350],[255 305 335],[240 310 345],[250 300 350],[245 300 345],[260 300 340]}, ...
 'macro_Q',{[45 65 95 130],[110 75 58 48],[55 95 58 110],[75 62 85 105],[36 48 62 80],[140 115 95 85]}, ...
 'carbonate',{false,true,false,true,false,true}, ...
 'fault_trace',{60,66,55,63,58,62},'throw',{31,28,35,30,36,25});
% Scenes: Q0 at 14 and 8 dB; Q0 x 0.7 at 14 dB on the same noise as the Q0 14 dB scene.
P.scenes=struct('tag',{'Q100_14dB','Q100_08dB','Q070_14dB'},'snr_db',{14,8,14},'q_scale',{1,1,.7},'level',{1,2,1});
P.seed=@(m,level,r) 1400000000+10000000*m+100000*level+r;
end
