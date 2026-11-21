"""Toy, field convention: is an overestimated SNR enough to push compensation below R0 on the box-band tie?
Same model as t_convention.py (imported by path). Operator noise power = true noise / kappa, kappa in {1, 10, 100}.
Noise is set so the scoring window is noise-limited at the upper band (shallow SNR 6 dB and 0 dB)."""
import importlib.util, sys, json
import numpy as np
spec = importlib.util.spec_from_file_location('tc', sys.argv[1]); tc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tc)
out = {}
for snr in [6.0, 0.0]:
    for kappa in [1.0, 10.0, 100.0]:
        rng = np.random.default_rng(11)
        acc = {k: {m: [] for m in tc.MS} for k in ['R0', 'UQ25', 'AUQ', 'UQ100']}
        for rep in range(6):
            Dn, ref, npow, spow = tc.synth(rng, 25.0, snr, False)
            S0 = spow * tc.W ** 2; Sa = S0 * np.exp(-2 * np.pi * tc.fs * tc.TA / tc.Q)
            y0 = np.fft.irfft(Dn * tc.band, tc.N); box = tc.r0_band(y0)
            for m in tc.MS:
                ys = {'R0': y0,
                      'UQ25': tc.apply(Dn, tc.t / tc.Q, 25.0, m, S0, npow / kappa),
                      'UQ100': tc.apply(Dn, tc.t / tc.Q, 100.0, m, S0, npow / kappa),
                      'AUQ': tc.apply(Dn, np.maximum(tc.t - tc.TA, 0) / tc.Q, 25.0, m, Sa, npow / kappa)}
                for k, y in ys.items():
                    acc[k][m].append(tc.score(y, ref, box)[2])
        key = f'snr{int(snr)}_kappa{int(kappa)}'
        out[key] = {k: {str(m): float(np.mean(v)) for m, v in d.items()} for k, d in acc.items()}
        print(f'--- shallow SNR {snr:g} dB, operator noise = true/{kappa:g}   TIEbox   m = ' + ' '.join(f'{m:g}' for m in tc.MS))
        for k in ['R0', 'UQ100', 'UQ25', 'AUQ']:
            print(f'   {k:6s} ' + '  '.join(f'{out[key][k][str(m)]:.3f}' for m in tc.MS))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
