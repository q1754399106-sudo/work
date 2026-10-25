from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];src=(ROOT.parent/'complex_fieldlike_20261023/core_p1/run_complex_model.m').read_text(encoding='utf-8')
src=src.replace('function run_complex_model(outdir)','function run_denoise_trial(outdir,datasetName)')
src=src.replace("'complex_input.json'","'denoised_input.json'")
src=src.replace("p.t=p.t(:);assert(~any(isfield(p,{'Q','fpeak','truth','rho','vs'})));","p.t=p.t(:);assert(~any(isfield(p,{'Q','fpeak','truth','rho','vs'})));\ninput.datasets=input.datasets(strcmp({input.datasets.name},datasetName));assert(numel(input.datasets)==1);")
start=src.index('policy=jsondecode(');end=src.index('offsets=repmat',start)
src=src[:start]+"saved=jsondecode(fileread(fullfile(root,'inputs','shared_calibration.json')));\nest=saved.estimate;corrected=saved.corrected;estbase=saved.original_estimate;\nfor field={'vp','h','offset','lambda_grid'},f=field{1};corrected.(f)=reshape(corrected.(f),1,[]);end\ncorrected.t=corrected.t(:);cal=saved.diagnostics;\ncmfull=uq_inverse_config(corrected,est);cmbase=uq_inverse_config(corrected,estbase);\n"+src[end:]
(ROOT/'core_p1/run_denoise_trial.m').write_text(src,encoding='utf-8')
for label,name in [('moderate','Moderate'),('strong','Strong')]:
 (ROOT/'core_p1'/('run_denoise_'+label+'.m')).write_text(f"function run_denoise_{label}(outdir)\nrun_denoise_trial(outdir,'{name}');\nend\n",encoding='utf-8')
