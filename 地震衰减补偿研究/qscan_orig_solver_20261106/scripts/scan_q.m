function S = scan_q(d, L, fw)
% Scan constant Q with the fixed inversion wavelet fw (first L.nscan traces). On every grid point:
%  C1 eq. 17 of Cheng et al. 2025 (CJG), signed: (sum|D_shallow| - sum|D_deep|) / sum|D_shallow| of the reconstructed data in
%     the band, trace-average spectra of the first and the last Tcrit seconds. Selection: scanning from weak compensation, the
%     first minimum of |eq. 17|: the + to - crossing (linear in v) or, without a crossing, the vertex of the parabola through the
%     first local minimum of |E1|; the unstable region beyond it is not used.
%  C2 Cheng et al. 2019 (SEG): l_0.1 norm of the inverted reflectivity, each trace scaled by its maximum. Selection: the most
%     prominent interior local minimum, vertex of the parabola in v (the global minimum falls on the strong-compensation edge).
d=d(:,1:L.nscan); nv=numel(L.vgrid); E1=nan(1,nv); C2=nan(1,nv); nT=round(L.Tcrit/L.dt);
nf=round(1/L.dt); f=(0:nf-1)'/(nf*L.dt); bm=f>=L.band(1)&f<=L.band(2);
for k=1:nv
    [r,comp]=qcdn_core(d,L.dt,fw,L.liter,L.niter,1/L.vgrid(k));
    Ds=mean(abs(fft(comp(1:nT,:),nf)),2); Dd=mean(abs(fft(comp(end-nT+1:end,:),nf)),2);
    E1(k)=(sum(Ds(bm))-sum(Dd(bm)))/sum(Ds(bm));
    rn=r./max(abs(r),[],1); C2(k)=mean(sum(abs(rn).^0.1,1));
end
S.E1=E1; S.C2=C2; S.v1=crossing(L.vgrid,E1); S.v2=prominent(L.vgrid,C2);
end
function v=crossing(g,e)
v=NaN;
for k=1:numel(e)-1
    if e(k)>0 && e(k+1)<=0
        v=g(k)+(g(k+1)-g(k))*e(k)/(e(k)-e(k+1)); return
    end
    if k>1 && abs(e(k))<abs(e(k-1)) && abs(e(k+1))>=abs(e(k))
        v=g(k); a=polyfit(g(k-1:k+1),abs(e(k-1:k+1)),2);
        if a(1)>0, v=min(max(-a(2)/(2*a(1)),g(k-1)),g(k+1)); end
        return
    end
end
end
function v=prominent(g,J)
best=-inf; v=NaN;
for m=2:numel(g)-1
    if J(m)<J(m-1) && J(m)<=J(m+1)
        p=min(max(J(1:m-1)),max(J(m+1:end)))-J(m);
        if p>best
            best=p; a=polyfit(g(m-1:m+1),J(m-1:m+1),2); v=g(m);
            if a(1)>0, v=min(max(-a(2)/(2*a(1)),g(m-1)),g(m+1)); end
        end
    end
end
end
