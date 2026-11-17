function [Tt,tt]=two_way(F,k,P)
% Two-way traveltime and t* for shot k at every node: the source leg plus the leg to the specular
% receiver of a flat reflector, x_r = 2x - x_s, taken from the shot-position fields by linear
% interpolation in x_r (clamped to the shot range).
nz=P.nz;nx=P.nx;ns=numel(P.ixs);xs=P.x(P.ixs);
q=(2*P.x-xs(k)-xs(1))/(xs(2)-xs(1))+1;q=min(max(q,1),ns);k0=floor(q);k0(k0==ns)=ns-1;w=q-k0;
base=(0:nx-1)*nz;i=(1:nz)';
i0=i+base+(k0-1)*nz*nx;i1=i0+nz*nx;w=repmat(w,nz,1);
Tt=F.T(:,:,k)+(1-w).*F.T(i0)+w.*F.T(i1);
tt=F.ts(:,:,k)+(1-w).*F.ts(i0)+w.*F.ts(i1);
end
