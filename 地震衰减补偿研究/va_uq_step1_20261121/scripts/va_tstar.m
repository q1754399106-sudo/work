function TS=va_tstar(M,P,C,vb,lin,xcol)
% Two-way t* per image node and shot (specular receiver, two_way.m) in the band velocity vb with the
% true compensation Q (E1, scaled by 1/m later), and the anchored difference (PLAN_ZH.md section 4):
%   Delta t* = t* accumulated below the anchor along the same ray
%            = t* computed with 1/Q_c = 0 above z_a(x) (Fermat: dT/de at the stationary ray, so only
%              the part of the ray below the anchor contributes),
% z_a(x) = M.spec.iface{C.anchor_iface}(x) (macro interface). This equals t* - t*_above-anchor of the
% same source-receiver ray. (Subtracting the t* of the anchor node at the same x and the same
% source-receiver pair is not the same ray: at large offsets it crosses the shallow low-Q layer more
% obliquely; smoke check: D4 W_deep Delta/t* median 0.05 that way.) Fields cached as
% cache/fields_<name>_vb.mat and cache/fields_<name>_vb_anc.mat. xcol: grid column of each node in lin.
F=tstar_fields(M,P,M.Qm,vb,'_vb');ns=numel(P.ixs);
[X,Z]=meshgrid(P.x,P.z);za=M.spec.iface{C.anchor_iface}(X);Qa=M.Qm;Qa(max(Z,0)<za)=Inf;
Fa=tstar_fields(M,P,Qa,vb,'_vb_anc');
TS.abs=zeros(numel(lin),ns,'single');TS.dlt=TS.abs;TS.za=za(1,:);TS.xcol=xcol;
for s=1:ns
 [~,tt]=two_way(F,s,P);TS.abs(:,s)=tt(lin);
 [~,tt]=two_way(Fa,s,P);TS.dlt(:,s)=max(tt(lin),0);
end
end
