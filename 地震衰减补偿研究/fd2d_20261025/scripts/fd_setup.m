function P=fd_setup()
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
P.models=struct('name',{'H1','H2','H3','H4','H5','H6'}, ...
 'seed',{2026102701,2026102711,2026102721,2026102731,2026102741,2026102751}, ...
 'macro_vp',{[1840 2290 2760 3310],[1880 2330 2790 3340],[1790 2240 2810 3370],[1820 2360 2730 3270],[1860 2270 2780 3300],[1800 2320 2770 3280]}, ...
 'macro_h',{[245 295 355],[235 315 335],[255 300 350],[240 310 340],[250 290 360],[260 305 345]}, ...
 'macro_Q',{[55 75 100 135],[95 65 50 45],[45 85 55 115],[65 55 75 95],[35 50 65 85],[125 105 85 75]}, ...
 'carbonate',{false,true,false,true,false,true}, ...
 'fault_trace',{60,66,55,70,62,58},'throw',{25,35,30,40,28,32});
% Scenes: Q0 at 14 and 8 dB; Q0 x 0.7 at 14 dB on the same noise as the Q0 14 dB scene.
P.scenes=struct('tag',{'Q100_14dB','Q100_08dB','Q070_14dB'},'snr_db',{14,8,14},'q_scale',{1,1,.7},'level',{1,2,1});
P.seed=@(m,level,r) 1400000000+10000000*m+100000*level+r;
end
