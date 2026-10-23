function build_demo_inputs(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));cache=plan.cache;
z=load(fullfile(cache,'cases','OriginalReference_WeakSection_12dB_20261411.mat'));cal=load(fullfile(cache,'calibration_input_OriginalReference.mat'));
D=z.D;offsets=z.offsets;public=cal.public;calibrationD=cal.calD;calibrationOffsets=cal.public.offset;domain='Section';
save(fullfile(outdir,'example_inputs.mat'),'D','offsets','public','calibrationD','calibrationOffsets','domain','-v7');
out=compensate_seismic(D,offsets,public,calibrationD,calibrationOffsets,domain);
save(fullfile(outdir,'example_output.mat'),'out','-v7');
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'optimization_flags',out.optimization_flags,'normal_residual',out.normal_residual,'input_fields',{{'D','offsets','public','calibrationD','calibrationOffsets','domain'}}));
end
