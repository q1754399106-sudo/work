function test_joint_api(outdir)
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);
[configs,~]=p1_models('development');c=configs{1};public=uq_public_config(c,c.nominal_frequency);
models=cell(1,11);for j=1:11,[a,b]=ray_absorption(c.t,c.offset(j),c);models{j}=struct('G',attenuation_matrix(a+b,c));end
calScene=make_stability_scene('SmoothAVO',c);[calD,~]=uq_observation(calScene,c,models,1:11,20261201,18);
scene=make_stability_scene('WeakSection',c);[D,~]=uq_observation(scene,c,models,6*ones(1,11),20261303,18);
offsets=repmat(c.offset(6),1,11);options=struct('mode','Physical','beta_factor',.3);
out=joint_compensate(D,offsets,public,calD,c.offset,'Section',options);
a=uq_audit_outputs(out);assert(a.flags==0&&a.max_normal<2e-7);assert(all(isfinite(out.gather(:))));
bad=public;bad.Q=c.Q;rejected=false;
try,joint_compensate(D,offsets,bad,calD,c.offset,'Section',options);catch ex,rejected=contains(ex.message,'Nonpublic');end
assert(rejected,'Truth-bearing public config was not rejected');
save(fullfile(outdir,'api_input_output.mat'),'D','public','calD','offsets','options','out','-v7');
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'finite_output',true,'truth_field_rejected',rejected, ...
 'solver_flags',a.flags,'max_normal_residual',a.max_normal,'shared_regularizer_difference',out.regularizer_difference));
fprintf('Research API end-to-end check passed; normal %.3e\n',a.max_normal);
end
function restore(p,r),path(p);rng(r);end
