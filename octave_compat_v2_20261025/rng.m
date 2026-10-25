function rng(seed,generator)
% Octave shim (v2): rng(seed,'twister') starts a MATLAB-compatible mt19937ar stream shared by
% rand.m and randn.m. Words come from mlwords.py and are buffered in Octave.
global ML2
assert(nargin==2&&strcmp(generator,'twister')&&isscalar(seed)&&seed>=0&&seed==fix(seed),'rng shim supports rng(seed,''twister'') only.');
ML2=struct('seed',double(seed),'W',[],'k',1);ml2_words(2^20);
end
