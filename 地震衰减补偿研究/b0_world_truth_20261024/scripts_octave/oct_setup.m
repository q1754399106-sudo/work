function P=oct_setup()
% Octave path: compatibility shims first (rng/randn = MATLAB R2024b twister + Ziggurat, sum 'all'
% and vecdim, string), then the unchanged scripts/ of this study, then bw_setup as in MATLAB.
here=fileparts(mfilename('fullpath'));root=fileparts(here);
addpath(fullfile(fileparts(root),'octave_compat_20261025'),'-begin');addpath(fullfile(root,'scripts'));
warning('off','all');P=bw_setup();
end
