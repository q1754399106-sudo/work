function o=structure_options()
root=fileparts(fileparts(mfilename('fullpath')));
o=jsondecode(fileread(fullfile(root,'inputs','options.json')));
end
