function F=tstar_fields(M,P,Qc,vmig,tag)
% Traveltime and t* from every shot position in the migration velocity and a compensation Q model.
% t* = dT/de of the traveltime in the slowness s (1 + e/Q) at e = 0 (Fermat: the ray is stationary),
% by the central difference e = +-0.5. Defaults (stage 1): Qc = M.Qm, vmig = M.vmig, tag ''.
% Cache: cache/fields_<name><tag>.mat.
if nargin<3||isempty(Qc),Qc=M.Qm;end
if nargin<4||isempty(vmig),vmig=M.vmig;end
if nargin<5,tag='';end
fn=fullfile(P.cache,['fields_' M.name tag '.mat']);
if isfile(fn),F=load(fn);return;end
s=1./vmig;e=.5;ns=numel(P.ixs);F.T=zeros(P.nz,P.nx,ns,'single');F.ts=F.T;
for k=1:ns
 F.T(:,:,k)=eikonal_fs(s,P.dx,P.izs,P.ixs(k));
 F.ts(:,:,k)=(eikonal_fs(s.*(1+e./Qc),P.dx,P.izs,P.ixs(k))-eikonal_fs(s.*(1-e./Qc),P.dx,P.izs,P.ixs(k)))/(2*e);
 if mod(k,10)==0,fprintf('fields %s%s shot %d/%d\n',M.name,tag,k,ns);end
end
save(fn,'-struct','F','-v7.3');
end
