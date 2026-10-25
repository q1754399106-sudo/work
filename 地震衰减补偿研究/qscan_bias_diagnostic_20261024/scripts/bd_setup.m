function P=bd_setup()
% Paths and constants. Every input is read from cache/inputs, a byte-verified copy of Git
% objects (scripts/snapshot_inputs.py), so shared sparse-checkout changes cannot break a run.
P.root=fileparts(fileparts(mfilename('fullpath')));P.in=fullfile(P.root,'cache','inputs');
assert(isfile(fullfile(P.root,'cache','SNAPSHOT.json')),'Run scripts/snapshot_inputs.py first.');
bae=fullfile(P.in,'bae_q_20261024');addpath(fullfile(bae,'core'),fullfile(bae,'source'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid, verbatim
P.K=6;P.threads=1; P.window=[.02 .25];        % bootstrap draws; noise window of estimate_trace_snr
end
