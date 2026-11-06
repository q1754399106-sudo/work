function S = scan_q2(d, L, fw)
% scan_q.m of qscan_orig_solver_20261106 (same solver, same eq. 17 and C2 rules), adding chengD on the same solutions:
% |slope| over 10-60 Hz of ln(A_deep/A_shallow), A = RMS amplitude spectrum over traces of the first and last third of the record
% (Hann taper, nfft 256), as crit_kappa.py of qscan_anchor_theory_20261030; grid minimum, then parabola vertex in v.
d=d(:,1:L.nscan); nv=numel(L.vgrid); E1=nan(1,nv); C2=nan(1,nv); CD=nan(1,nv); nT=round(L.Tcrit/L.dt);
nf=round(1/L.dt); f=(0:nf-1)'/(nf*L.dt); bm=f>=L.band(1)&f<=L.band(2);
for k=1:nv
    [r,comp]=qcdn_core(d,L.dt,fw,L.liter,L.niter,1/L.vgrid(k));
    Ds=mean(abs(fft(comp(1:nT,:),nf)),2); Dd=mean(abs(fft(comp(end-nT+1:end,:),nf)),2);
    E1(k)=(sum(Ds(bm))-sum(Dd(bm)))/sum(Ds(bm));
    rn=r./max(abs(r),[],1); C2(k)=mean(sum(abs(rn).^0.1,1));
    CD(k)=cheng_d(comp,L);
end
S.E1=E1; S.C2=C2; S.CD=CD; [S.v1,S.cross]=crossing(L.vgrid,E1); S.vd=gridmin(L.vgrid,CD); S.v2=prominent(L.vgrid,C2);
end
function s=cheng_d(comp,L)
n=size(comp,1); m=round(n/3); h=0.5-0.5*cos(2*pi*(0:m-1)'/(m-1)); nf=256; f=(0:nf/2)'/(nf*L.dt);
As=fft(comp(1:m,:).*h,nf); Ad=fft(comp(end-m+1:end,:).*h,nf);
As=sqrt(mean(abs(As(1:nf/2+1,:)).^2,2)); Ad=sqrt(mean(abs(Ad(1:nf/2+1,:)).^2,2)); b=f>=10&f<=60;
a=polyfit(f(b),log(max(Ad(b),1e-30)./max(As(b),1e-30)),1); s=abs(a(1));
end
function v=gridmin(g,c)
[~,m]=min(c); v=g(m);
if m>1&&m<numel(g)
    a=polyfit(g(m-1:m+1),c(m-1:m+1),2);
    if a(1)>0, v=min(max(-a(2)/(2*a(1)),g(m-1)),g(m+1)); end
end
end
function [v,cross]=crossing(g,e)
v=NaN; cross=false;
for k=1:numel(e)-1
    if e(k)>0 && e(k+1)<=0
        v=g(k)+(g(k+1)-g(k))*e(k)/(e(k)-e(k+1)); cross=true; return
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
