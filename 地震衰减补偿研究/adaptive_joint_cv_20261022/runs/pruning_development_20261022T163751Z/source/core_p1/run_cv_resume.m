function run_cv_resume(outdir)
root=fileparts(mfilename('fullpath'));a=jsondecode(fileread(fullfile(root,'inputs','cv_resume.json')));
run_cv_screen(outdir,9,true,a.path);
end
