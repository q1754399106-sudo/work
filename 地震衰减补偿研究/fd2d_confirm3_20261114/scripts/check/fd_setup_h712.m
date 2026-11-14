function P=fd_setup_h712()
% fd2d_confirm_20261109: fd2d_20261025/scripts/fd_setup.m with new models H7-H12 only (UQ-W confirmation).
% Paths, constants, the six 2-D development geologies H1-H6 and the three scenes of PLAN_ZH.md.
% Every upstream input is read from cache/inputs, a byte-verified copy of Git objects (snapshot_inputs.py).
P.root=fileparts(fileparts(fileparts(mfilename('fullpath'))));P.in=fullfile(P.root,'cache','inputs');
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
P.models=struct('name',{'H7','H8','H9','H10','H11','H12'}, ...
 'seed',{2026110901,2026110911,2026110921,2026110931,2026110941,2026110951}, ...
 'macro_vp',{[1830 2300 2750 3290],[1870 2320 2800 3350],[1810 2260 2790 3360],[1850 2340 2740 3280],[1840 2280 2770 3310],[1820 2310 2780 3290]}, ...
 'macro_h',{[250 300 345],[240 310 340],[245 305 355],[255 295 345],[235 300 350],[250 310 340]}, ...
 'macro_Q',{[50 70 105 140],[100 70 55 45],[50 90 60 120],[70 60 80 100],[38 52 68 88],[130 110 90 80]}, ...
 'carbonate',{false,true,false,true,false,true}, ...
 'fault_trace',{64,57,68,53,61,65},'throw',{27,33,29,38,26,34});
% Scenes: Q0 at 14 and 8 dB; Q0 x 0.7 at 14 dB on the same noise as the Q0 14 dB scene.
P.scenes=struct('tag',{'Q100_14dB','Q100_08dB','Q070_14dB'},'snr_db',{14,8,14},'q_scale',{1,1,.7},'level',{1,2,1});
P.seed=@(m,level,r) 1400000000+10000000*m+100000*level+r;
end
