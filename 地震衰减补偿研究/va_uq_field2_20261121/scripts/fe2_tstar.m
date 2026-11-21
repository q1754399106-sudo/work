function TS=fe2_tstar(C,P,lin)
% Two-way t* per interior node and shot (nominal Q_m = 120; scaled by 1/m when the operators are formed).
%   TS.abs : absolute t*, the stage-4 cache uq_rtm_20261117/cache/fld/ts.mat written by fld_tstar.m (reused,
%            read only; the operators UQW1_100, UQW1_25 use it as in stage 4 and its post hoc check).
%   TS.dlt : anchored Delta t* (VA-UQ; same-ray definition of va_uq_step1_20261121): fld_tstar.m unchanged,
%            with 1/Q = 0 above the anchor z_a(x) (cache/anchor.mat), its own cache in cache/anc/.
%            Fermat: dT/de is taken on the stationary ray, so only the ray length below z_a contributes;
%            Delta t* = t* - t*(part of the same ray above the anchor). Clipped at 0 (rounding).
S=load(fullfile(C.s4cache,'ts.mat'));TS.abs=S.TS;clear S
A=load(fullfile(C.cache,'anchor.mat'));Pa=P;Pa.cache=fullfile(C.cache,'anc');if ~isfolder(Pa.cache),mkdir(Pa.cache);end
Z=repmat(P.z(:),1,P.nx);Qa=P.Qm;Qa(Z<repmat(A.za_full,P.nz,1))=Inf;Pa.Qm=Qa;
TS.dlt=max(fld_tstar(Pa,lin),0);TS.za_full=A.za_full;
end
