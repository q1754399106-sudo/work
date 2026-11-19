function TS=fld_tstar(P,lin)
% Two-way t* (nominal Q model) for each shot at every interior node: eikonal dT/de at every shot position
% (as tstar_fields.m), source leg plus the leg to the flat-reflector specular receiver x_r = 2x - x_s,
% interpolated between shot positions (as two_way.m). Cached in cache/fld/ts.mat.
fn=fullfile(P.cache,'ts.mat');if isfile(fn),S=load(fn);TS=S.TS;return;end
s=1./P.vmig;e=.5;ns=P.nshot;F=zeros(P.nz,P.nx,ns,'single');
for k=1:ns
 F(:,:,k)=(eikonal_fs(s.*(1+e./P.Qm),P.dx,P.src_iz(k),P.src_ix(k))-eikonal_fs(s.*(1-e./P.Qm),P.dx,P.src_iz(k),P.src_ix(k)))/(2*e);
 if mod(k,20)==0,fprintf('t* fields %d/%d\n',k,ns);end
end
[xs,o]=sort(P.src_x);TS=zeros(numel(lin),ns,'single');[iz,ix]=ind2sub([P.nz P.nx],lin);
for k=1:ns
 xr=2*P.x(ix)'-P.src_x(k);q=interp1(xs,1:ns,min(max(xr,xs(1)),xs(end)));k0=min(floor(q),ns-1);w=q-k0;
 a=sub2ind([P.nz P.nx ns],iz,ix,o(k0));b=sub2ind([P.nz P.nx ns],iz,ix,o(k0+1));
 TS(:,k)=F(sub2ind([P.nz P.nx ns],iz,ix,k*ones(size(iz))))+single(1-w).*F(a)+single(w).*F(b);
end
save(fn,'TS','-v7.3');
end
