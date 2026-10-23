function profile_one(root,outdir,preflight)
addpath(fullfile(root,'core_p1','matlab_src')); mkdir(outdir);
[cc,~]=p1_models('historical');c=cc{1};p=uq_public_config(c,c.nominal_frequency);
est=jsondecode(fileread(fullfile(preflight,'SoftLowQ','matlab','calibration_20261030.json')));
cm=uq_inverse_config(p,est);models=prepare_models_v2(c.offset,cm);R=avo_projector_v2(models);
trueModels=cell(1,11);
for j=1:11, [a,b]=ray_absorption(c.t,c.offset(j),c); trueModels{j}=struct('G',attenuation_matrix(a+b,c)); end
[mm,domain,map,scene]=uq_case('SmoothAVO',c,models);
[D,~]=uq_observation(scene,c,trueModels,map,20261101,18);
names={'Baseline','Frozen','StrongerDamping','NoiseUncertaintyOnly','QModelErrorOnly','CombinedUncertainty','CombinedUncertainty_NoStructure'};
profile on
times=zeros(size(names));
for i=1:numel(names)
 tic;[o,~,~]=p1_solve_branch(D,mm,map,cm,domain,R,names{i});times(i)=toc;
 fprintf('%s %.3fs %.3e\n',names{i},times(i),o.direct_record.normal_residual);
end
profile off; info=profile('info');
[~,idx]=sort([info.FunctionTable.TotalTime],'descend');
top=info.FunctionTable(idx(1:min(40,numel(idx))));
uq_writejson(fullfile(outdir,'profile.json'),top);
uq_writejson(fullfile(outdir,'method_times.json'),struct('names',{names},'seconds',times));
end
