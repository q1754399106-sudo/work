from pathlib import Path
import zipfile,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];archive=ROOT.parent/'optimized_seismic_method_20261023.zip';dest=ROOT.parent/'delivery_check_20261023';assert not dest.exists();dest.mkdir()
with zipfile.ZipFile(archive) as z:
    for name in z.namelist():assert not name.startswith(('/','\\')) and '..' not in Path(name).parts
    z.extractall(dest)
q=lambda p:p.as_posix().replace("'","''")
core=dest/'core_p1';output=dest/'output.mat'
code=f"addpath('{q(core)}');a=load('{q(core/'examples/example_inputs.mat')}');out=compensate_seismic(a.D,a.offsets,a.public,a.calibrationD,a.calibrationOffsets,a.domain);assert(out.validation.ready_for_research_use&&out.optimization_flags==0);save('{q(output)}','out','-v7');disp('PORTABLE_PACKAGE_CHECK_PASSED');"
cmd=['D:/Program Files (x86)/MATLAB2024B/bin/matlab.exe','-wait','-nojvm','-singleCompThread','-batch',code];start=time.monotonic()
with (dest/'stdout.log').open('wb') as stdout,(dest/'stderr.log').open('wb') as stderr:
    p=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    try:exit_code=p.wait(timeout=600)
    except subprocess.TimeoutExpired:
        subprocess.run(['taskkill','/PID',str(p.pid),'/T','/F'],capture_output=True);raise
result={'passed':exit_code==0 and output.exists(),'exit_code':exit_code,'unpacked_directory':str(dest),'elapsed_seconds':time.monotonic()-start,'test':'Unpack code ZIP into a fresh directory and run the public entry using only bundled observation/calibration inputs; require ready research status and no optimization failures.'}
(ROOT/'analysis/delivery_check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False,indent=2));assert result['passed']
