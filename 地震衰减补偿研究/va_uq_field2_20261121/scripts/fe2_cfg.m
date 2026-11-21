function [C,P]=fe2_cfg()
% Configuration of va_uq_field2_20261121 (PLAN_ZH.md): VA-UQ, second field pass on the Volve line of
% uq_rtm_20261117 stage 4. Engine, data, geometry, evaluation are the stage-4 functions (addpath,
% unchanged); the VA-UQ table va_tab.m is the one of va_uq_step1_20261121 (addpath, unchanged).
% Stage-4 caches are read only; everything computed here goes to this experiment's cache/.
here=fileparts(mfilename('fullpath'));C.root=fileparts(here);base=fileparts(C.root);
C.s4=fullfile(base,'uq_rtm_20261117');C.st1=fullfile(base,'va_uq_step1_20261121');
for d={fullfile(C.s4,'scripts'),fullfile(C.st1,'scripts')}
 if ~contains(path,d{1}),addpath(d{1});end
end
P=fld_setup();                                   % stage-4 constants, data, velocity, Q_m = 120, W
C.s4cache=P.cache;                               % uq_rtm_20261117/cache/fld (read only)
C.cache=fullfile(C.root,'cache');if ~isfolder(C.cache),mkdir(C.cache);end
C.ms=[.5 .6 .7 .85 1 1.2 1.4 1.7 2];             % stage-4 m grid, Q_c = m * 120
C.sig=P.sig;                                     % 0.447
C.fv=[15 35];                                    % VA-UQ: ln f_v uniform on [15, 35] Hz
C.B=300;C.flip_seed=20261118;                    % block flip as est_block.m: rng(flip_seed + B)
C.anchor=struct('zwin',[2150 2400],'jump_half',20,'medlen',21);   % PLAN section 3 (fe2_anchor.m)
C.threads=12;
C.mem_need_gb=NaN;                               % set from the smoke test (fe2_run memory gate)
f=fullfile(C.cache,'mem_need_gb.txt');if isfile(f),C.mem_need_gb=str2double(fileread(f));end
% methods formed in pass 2: name -> {table, t* key, noise estimate}
% (VAfix, f_v known, was dropped after the migration smoke: 255 s per frequency for pass 2 on the loaded
% machine; every table family costs three interpolations, an extra method on an existing family only one H)
C.meth={'VA','dlt','va','B300';'UQW1_100','abs','u100','B300';'UQW1_25','abs','u25','B300'; ...
        'VA_B1','dlt','va','B1';'UQW1_100_B1','abs','u100','B1';'UQW1_25_B1','abs','u25','B1'};
end
