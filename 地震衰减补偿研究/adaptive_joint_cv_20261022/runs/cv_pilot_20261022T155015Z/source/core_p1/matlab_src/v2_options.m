function o=v2_options()
root=fileparts(fileparts(mfilename('fullpath')));
o=jsondecode(fileread(fullfile(root,'inputs','v2_options.json')));
end
