function out=compensate_seismic(D,offsets,public,calibrationD,calibrationOffsets,domain)
% Single user entry: the same FullCalRefit rule for every model.
out=reflectivity_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain);
root=fileparts(mfilename('fullpath'));p=fullfile(root,'inputs','validation_status.json');
out.research_only=true;
if isfile(p),out.validation=jsondecode(fileread(p));end
end
