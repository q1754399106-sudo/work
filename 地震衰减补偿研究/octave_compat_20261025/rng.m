function rng(seed,generator)
% Octave shim: rng(seed,'twister') starts a MATLAB-compatible mt19937ar stream for randn.m.
global ML_RNG_STATE
assert(nargin==2&&strcmp(generator,'twister')&&isscalar(seed)&&seed>=0&&seed==fix(seed),'rng shim supports rng(seed,''twister'') only.');
ML_RNG_STATE=struct('seed',double(seed),'offset',0);
end
