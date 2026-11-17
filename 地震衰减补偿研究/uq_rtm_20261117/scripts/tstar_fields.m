function F=tstar_fields(M,P)
% Traveltime and t* from every shot position in the migration velocity and the compensation Q model.
% t* = dT/de of the traveltime in the slowness s (1 + e/Q) at e = 0 (Fermat: the ray is stationary),
% by the central difference e = +-0.5.
fn=fullfile(P.cache,['fields_' M.name '.mat']);
if isfile(fn),F=load(fn);return;end
s=1./M.vmig;e=.5;ns=numel(P.ixs);F.T=zeros(P.nz,P.nx,ns,'single');F.ts=F.T;
for k=1:ns
 F.T(:,:,k)=eikonal_fs(s,P.dx,P.izs,P.ixs(k));
 F.ts(:,:,k)=(eikonal_fs(s.*(1+e./M.Qm),P.dx,P.izs,P.ixs(k))-eikonal_fs(s.*(1-e./M.Qm),P.dx,P.izs,P.ixs(k)))/(2*e);
 fprintf('fields %s shot %d/%d\n',M.name,k,ns);
end
save(fn,'-struct','F','-v7.3');
end
