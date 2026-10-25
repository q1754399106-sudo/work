"""Independent small dense numerical checks; NOT seismic effectiveness evidence."""
from pathlib import Path
import argparse
import json
import numpy as np


def main(output):
    output=Path(output)
    output.mkdir(parents=True, exist_ok=False)
    rng=np.random.default_rng(20261019)
    n=24
    q=rng.normal(size=(n,5))*.12
    diag=np.exp(rng.normal(size=n))*.04
    base=np.diag(diag)
    cov=base+q@q.T
    precision=np.linalg.solve(cov,np.eye(n))
    bi=np.diag(1/diag)
    woodbury=bi-bi@q@np.linalg.solve(np.eye(q.shape[1])+q.T@bi@q,q.T@bi)
    transform=np.eye(n)+.08*rng.normal(size=(n,n))
    g=rng.normal(size=(n,n))/np.sqrt(n)
    truth=rng.normal(size=n)
    bias=rng.normal(size=n)*.2
    data=g@truth+bias
    centered=data-bias
    propagated=transform@cov@transform.T
    ps=transform.T@np.linalg.solve(propagated,transform)
    penalty=np.eye(n)*.02
    direct=np.linalg.solve(g.T@precision@g+penalty,g.T@precision@centered)
    split=np.linalg.solve(g.T@ps@g+penalty,g.T@ps@centered)
    # The floor is numerical independent noise, not measured model error.
    floor=.003
    floored=propagated+floor*np.eye(n)
    whiten=np.linalg.inv(np.linalg.cholesky(floored))
    expectation=np.trace(whiten@propagated@whiten.T)
    eigen=np.linalg.eigvalsh(propagated)
    analytic=np.sum(eigen/(eigen+floor))
    checks={
        'woodbury_relative_error':float(np.linalg.norm(woodbury-precision)/np.linalg.norm(precision)),
        'exact_propagation_precision_relative_error':float(np.linalg.norm(ps-precision)/np.linalg.norm(precision)),
        'matched_route_solution_relative_error':float(np.linalg.norm(direct-split)/np.linalg.norm(direct)),
        'bias_centering_residual':float(np.linalg.norm(g@truth-data+bias)),
        'wrong_bias_sign_residual':float(np.linalg.norm(g@truth-data-bias)),
        'floored_whitened_expected_energy':float(expectation),
        'floored_target_trace_error':float(abs(expectation-analytic)),
        'minimum_covariance_eigenvalue':float(np.linalg.eigvalsh(cov).min()),
        'ce_zero_precision_error':float(np.linalg.norm(np.linalg.solve(base,np.eye(n))-bi)),
        'shared_q_cross_trace_covariance_norm':float(np.linalg.norm((q@q.T)[:8,8:16])),
    }
    assert checks['woodbury_relative_error']<1e-12
    assert checks['exact_propagation_precision_relative_error']<1e-12
    assert checks['matched_route_solution_relative_error']<1e-12
    assert checks['bias_centering_residual']<1e-12
    assert checks['wrong_bias_sign_residual']>1e-3
    assert checks['floored_target_trace_error']<1e-10
    assert 0<expectation<n
    assert checks['minimum_covariance_eigenvalue']>0
    assert checks['ce_zero_precision_error']<1e-12
    result={'status':'PASS','seed':20261019,'checks':checks,
            'scope':'Independent dense linear algebra in Python. No MATLAB code execution, Q calibration, geological test or performance estimate.'}
    (output/'checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    main(parser.parse_args().output)
