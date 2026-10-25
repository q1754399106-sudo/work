from pathlib import Path
import hashlib, json, shutil
from datetime import datetime, timezone

base = Path(__file__).resolve().parent
old = base / 'q_dense_gathers_20261023'
root = base / 'soft_refit_20261023'
root.mkdir(exist_ok=False)
for name in ['source', 'inputs']:
    shutil.copytree(old / name, root / name)
for name in ['results', 'analysis', 'figures', 'cache']:
    (root / name).mkdir()
files = {}
for folder in ['source', 'inputs']:
    for p in (root / folder).rglob('*'):
        if p.is_file():
            files[p.relative_to(root).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
(root / 'UPSTREAM_HASHES.json').write_text(json.dumps(files, indent=2), encoding='utf-8')
plan = dict(created_utc=datetime.now(timezone.utc).isoformat(),
    question='Can continuous group shrinkage preserve weak events better than hard BIC group deletion?',
    baseline=str(old), errors_percent=[-30,-10,0,10,30], traces=99,
    observations='Exact prior 99-trace observations, targets and frozen estimated sigma; no new synthesis.',
    methods=['FixedSplitSharedQ','FullCalRefit','PrunedFullCalRefit','SoftGroupRefit'],
    formula='For each existing 0.5/f support cluster: gain=E/(E+V), E=a_g^T Wgram_g a_g, V=sigma^2 tr(Wgram_g K_g K_g^T). K=(C^T C+ridge I)^-1 C^T. Wgram=(W/scale)^T(W/scale).',
    controls='Keep original dictionary, support discovery, ridge rule, noise estimates, Q/wavelet, continuous blocks of 11. No neighbor averaging, no hard group deletion, no truth at inference.',
    uncertainty='Conditional fixed-support white-noise propagation proxy; ignores selection uncertainty, colored-noise covariance and Q-model bias. Not a calibrated posterior.',
    selection='One frozen candidate; all old models are retrospective development benchmarks. Report gains and losses versus all three baselines, do not claim independent generalization.',
    acceptance='Preserve weak reflection better than hard pruning while inspecting overall error and amplification. If unsuccessful do not silently promote candidate.',
    runtime='MATLAB 4 threads, checkpoint each block; timeout 3600s with process-alive and log-growth monitoring.')
(root / 'PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
(root / 'RUN_LEDGER.json').write_text(json.dumps(dict(status='prepared',completed_blocks=0,expected_blocks=180,outputs_root=str(root)),indent=2),encoding='utf-8')
print(root)
