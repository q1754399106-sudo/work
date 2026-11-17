function P=rtm_setup()
% Constants of uq_rtm_20261117 (PLAN_STAGE1_ZH.md). Every array (model, wavefield, image) is the
% full grid including the absorbing strips; P.in selects the interior.
P.root=fileparts(fileparts(mfilename('fullpath')));P.cache=fullfile(P.root,'cache');
P.dx=4;P.npml=25;P.nzi=400;P.nxi=600;                       % interior 1600 m x 2400 m
P.nz=P.nzi+2*P.npml;P.nx=P.nxi+2*P.npml;
P.z=((1:P.nz)-P.npml-1)*P.dx;P.x=((1:P.nx)-P.npml-1)*P.dx;  % z = 0, x = 0 at the first interior node
P.in={P.npml+(1:P.nzi),P.npml+(1:P.nxi)};
P.izs=P.npml+3;                                             % sources and receivers at z = 8 m
P.ixr=P.npml+1+(0:2:P.nxi-1);                               % 300 receivers, 8 m, x = 0-2392 m
P.ixs=P.npml+1+(50:10:550);                                 % 51 shots, 40 m, x = 200-2200 m
% time-domain finite-difference modelling (fd2d_visco, GSLS, fref = 100 Hz as the UQ operators)
P.fd=struct('dt',4e-4,'nt',5000,'t0',.06,'fpeak',25,'fref',100);
P.dt=.004;P.nt=500;P.t=(1:P.nt)'*P.dt;                    % stored records: t = 4 ms ... 2.0 s
P.f=3:.5:55;                                                % migration frequencies (df = 1/2 s)
P.smooth_m=30;                                              % migration velocity: Gaussian, sigma 30 m
P.threads=12;
% UQ operators (as qunknown_comp_20261109: lognormal sigma 0.447, gain limits in dB)
P.sig=.447;P.m=[.5 .6 .7 .85 1 1.2 1.4 1.7 2];
P.snr_win=.4;P.snr_hop=.1;P.snr_min=1e-3;
P.noise_db=[-12 -30 -40];                                   % trace-level noise levels (prep_data; round 2)
P.snr_box=400;                                              % image-domain SNR window (m), migrate
P.routeC=struct('data','m30','m',[.7 1 1.4]);              % route C (data-domain UQ-W, gamma 0; round 2)
P.mig_sets={'att',[.7 1 1.4];'m12',[.7 1 1.4];'m30',[.5 .6 .7 .85 1 1.2 1.4 1.7 2];'m40',[.5 .7 1 1.4 2]};
P.gr=struct('base','m30','tag','g30','v',[700 900],'f0',10,'db',20);  % stage 2 coherent noise (prep_data)
P.snr_band=5;                                              % +-Hz band of the noise estimate (round 3)
P.flip_seed=20261118;                                       % random-polarity noise estimate (round 2)
end
