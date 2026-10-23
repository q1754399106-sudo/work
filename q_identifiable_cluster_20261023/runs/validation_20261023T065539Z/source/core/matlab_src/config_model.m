function c = config_model()
% All units: seconds, metres, Hz. No optional MATLAB toolbox required.
c.dt = 0.002; c.nt = 601; c.t = (0:c.nt-1)'*c.dt;
c.offset = 0:150:1500; % set 0:50:1500 for thesis-style dense sampling
c.vp = [2100 2250 2450 2800];
c.vs = [550 700 850 1000]; c.rho = [1.8 2.1 2.3 2.4];
c.Q = [40 50 70 100];
c.event_t = [0.38 0.65 0.90]; % chosen, NOT recovered thesis parameters
c.h = diff([0 c.event_t]).*c.vp(1:3)/2;
c.fpeak = 60; c.fref = 100; c.snr_db = 18; c.seed = 20261017;
c.nfft = 2^nextpow2(4*c.nt); % padded linear-domain crop
c.lambda_x = 1e-3; c.lambda_fixed = 1e-3;
c.lambda_grid = logspace(-7,0,41);
c.cov_floor = 1e-6; % relative to sigma^2
c.q_scales = [0.8 1 1.2];
c.output = fullfile(fileparts(mfilename('fullpath')),'results');
end
