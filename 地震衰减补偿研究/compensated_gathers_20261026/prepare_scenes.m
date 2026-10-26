function prepare_scenes()
% New four-layer display gathers; archived generator and forward operator.
root=fileparts(mfilename('fullpath'));pkg=fullfile(root,'source','compensation_package_20261024');
addpath(pkg,fullfile(pkg,'source'),fullfile(pkg,'core'));maxNumCompThreads(4);
out=fullfile(root,'cache','cases');truthdir=fullfile(root,'cache','truth');
if ~isfolder(out),mkdir(out);end;if ~isfolder(truthdir),mkdir(truthdir);end
diary(fullfile(root,'prepare_MATLAB.log'));cl=onCleanup(@()diary('off'));
c0=config_model();spec=jsondecode(fileread(fullfile(pkg,'inputs','p1_experiment_spec.json')));
c0.model_name='OriginalControl';configs={c0};splits={'historical','development','heldout_new_geology'};
for k=1:numel(splits)
 s=spec.models.(splits{k});for j=1:numel(s)
  c=c0;fields={'vp','vs','rho','Q','event_t','fpeak'};for a=1:numel(fields),ff=fields{a};c.(ff)=s(j).(ff)(:)';end
  c.h=diff([0 c.event_t]).*c.vp(1:3)/2;c.model_name=s(j).name;configs{end+1}=c; %#ok<AGROW>
 end
end
for k=1:numel(configs)
 c=configs{k};n=c.model_name;dest=fullfile(out,[n '.mat']);if isfile(dest),continue;end
 assert(max(c.event_t)<max(c.t),'Event outside time range');
 scene=make_stability_scene('FaultSection',c);target=scene.truth;
 [a,b]=ray_absorption(c.t,750,c);G=attenuation_matrix(a+b,c);clean=G*target;
 seed=2026102600+k;rng(seed,'twister');snr_db=14;sigma0=norm(clean,'fro')/sqrt(numel(clean))/10^(snr_db/20);
 D=clean+sigma0*randn(size(clean));sigma=sigma0*ones(1,size(D,2));config=c;config.solver_offset=750;
 save(dest,'config','D','sigma','seed','snr_db','-v7');save(fullfile(truthdir,[n '.mat']),'target','clean','-v7');
 fprintf('Prepared %s (%d x %d)\n',n,size(D));
end
% V3 has calibrated, uncertain Q: apply existing QScan rule without truth.
dest=fullfile(out,'V3_qscan.mat');
if ~isfile(dest)
 S=load(fullfile(root,'cache','inputs','coherence_soft_20261023','inputs','validation','V3.mat'));
 D=S.obs(1).D;sigma=S.obs(1).sigma;config=S.config;
 started=tic;outcome=q_compensate(D,config,'mode','qscan','sigma',sigma);
 factor=outcome.qscan_factor;config.Q=outcome.applied_Q;
 save(dest,'config','factor','outcome','-v7');fprintf('V3 QScan complete factor %.6f, %.1fs\n',factor,toc(started));
end
fprintf('PREPARE COMPLETE\n');
end
