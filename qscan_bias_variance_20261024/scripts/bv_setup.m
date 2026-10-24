function P=bv_setup()
% Paths, frozen upstream solver code and the four F1/F2 scenes of bae_q_20261024.
% Replicate seeds are disjoint from the archived noise seeds (2026102452/53/62/63).
% Code and F1/F2 data come from a hash-verified byte copy in cache/ when present
% (scripts/snapshot_upstream.py), otherwise from the repository folder itself.
P.root=fileparts(fileparts(mfilename('fullpath')));P.repo=fileparts(P.root);
P.bae=fullfile(P.root,'cache','bae_q_20261024');
if ~isfile(fullfile(P.bae,'SNAPSHOT.json')),P.bae=fullfile(P.repo,'bae_q_20261024');end
addpath(fullfile(P.bae,'core'),fullfile(P.bae,'source'));
P.grid=exp(linspace(log(.6),log(2.2),13));   % QScan grid, verbatim
P.factors=[1 P.grid];                         % slot 1 is FullCalRefit at the public Q0
P.replicates=100;P.threads=4;
P.scenes=struct('model',{'F1','F1','F2','F2'},'condition',{'Moderate','Strong','Moderate','Strong'}, ...
 'snr_db',{14,8,14,8},'archived_seed',{2026102452,2026102453,2026102462,2026102463});
for s=1:numel(P.scenes)
 P.scenes(s).seed_base=1e9+1e6*str2double(P.scenes(s).model(2))+1e5*(1+strcmp(P.scenes(s).condition,'Strong'));
end
end
