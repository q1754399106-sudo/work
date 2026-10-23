function out=compensate_seismic(D,offsets,public,calibrationD,calibrationOffsets,domain)
% Single user entry: the same PrunedFullCalRefit rule for every model.
out=pruned_compensate(D,offsets,public,calibrationD,calibrationOffsets,domain);
root=fileparts(mfilename('fullpath'));p=fullfile(root,'inputs','validation_status.json');
out.research_only=true;
if isfile(p),out.validation=jsondecode(fileread(p));out.development_only=~out.validation.ready_for_research_use;end
end
