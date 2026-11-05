function fw = fw_anchor(fp, ts, va)
% Ricker dominant frequency of the source implied by the measured dominant frequency fp at time ts under an assumed
% attenuation va = 1/Q_a, a(ts) = ts: 2/fp - 2 fp/fw^2 - pi ts va = 0 (Zhang and Ulrych 2002). va = 0: fw = fp (Cheng et al. 2025).
fw = 1/sqrt(1/fp^2 - pi*ts*va/(2*fp));
end
