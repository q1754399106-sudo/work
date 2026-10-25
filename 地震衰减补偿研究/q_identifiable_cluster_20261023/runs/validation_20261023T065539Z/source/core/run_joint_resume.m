function run_joint_resume(outdir)
root=fileparts(mfilename('fullpath'));
resume=jsondecode(fileread(fullfile(root,'inputs','resume_pointer.json')));
run_joint_trial(outdir,resume);
end
