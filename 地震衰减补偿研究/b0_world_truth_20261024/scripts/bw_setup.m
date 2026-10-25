function P=bw_setup()
% Paths and constants. Every input is read from cache/inputs, a byte-verified copy of Git
% objects (scripts/snapshot_inputs.py), so shared sparse-checkout changes cannot break a run.
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid, verbatim
P.K=6;P.threads=1;                            % as in qscan_bias_diagnostic_20261024 (B0 reproduces exactly)
P.coherence=struct('half_width',2,'dip_periods',.25,'lambda_min',.5,'lambda_max',16);   % round-3 "Standard"
end
