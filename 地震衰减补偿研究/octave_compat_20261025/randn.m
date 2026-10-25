function X=randn(varargin)
% Octave shim: MATLAB R2024b randn after rng(seed,'twister') (Ziggurat), via mlrandn.py.
global ML_RNG_STATE
assert(~isempty(ML_RNG_STATE),'randn shim: call rng(seed,''twister'') first.');
if nargin==0,dims=[1 1];elseif nargin==1,dims=double(varargin{1});if isscalar(dims),dims=[dims dims];end
else,dims=cellfun(@double,varargin);end
n=prod(dims);here=fileparts(mfilename('fullpath'));out=[tempname() '.bin'];
cmd=sprintf('python3 -I "%s" %d %d %d "%s"',fullfile(here,'mlrandn.py'),ML_RNG_STATE.seed,ML_RNG_STATE.offset,n,out);
[status,txt]=system(cmd);assert(status==0,'mlrandn.py failed: %s',txt);
fid=fopen(out,'r');X=fread(fid,n,'double',0,'ieee-le');fclose(fid);delete(out);
ML_RNG_STATE.offset=str2double(strtrim(txt));X=reshape(X,dims);
end
