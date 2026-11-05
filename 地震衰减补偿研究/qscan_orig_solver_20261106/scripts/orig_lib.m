function L = orig_lib()
% Shared definitions of the experiment (PLAN_ZH.md).
L.dt=0.002; L.n1=1001; L.f0=30; L.liter=10; L.niter=10;             % values of the original synthetic (qcomp/synth/SConstruct)
L.ntr=384; L.nscan=24;                                              % traces for the dominant frequency / traces inverted
L.t=(0:L.n1-1)'*L.dt;
L.vgrid=exp(linspace(log(1/400),log(1/20),21));                     % scanned v = 1/Q
L.Tcrit=0.25; L.band=[5 90];                                        % eq. 17 of Cheng et al. 2025 (CJG): window length and band
L.windows=struct('W1',[0 0.2],'W2',[0.4 0.6],'W3',[0.9 1.1]);
L.noise=0.05;
end
