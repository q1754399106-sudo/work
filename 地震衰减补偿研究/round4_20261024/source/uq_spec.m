function o=uq_spec()
root=fileparts(fileparts(mfilename('fullpath'))); o=jsondecode(fileread(fullfile(root,'inputs','uq_spec.json')));
end
