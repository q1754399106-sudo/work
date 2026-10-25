% Verbatim from qcriteria_20261024/run_criteria.m (function spectral_factor).
function c=spectral_factor(D,cm,offset,grid)
% Spectral-ratio t* differences between time windows vs dictionary prediction.
win=[.30 .60;.60 .90;.90 1.20;1.20 1.60];t=cm.t;N=512;f=(0:N/2)'/(N*cm.dt);band=f>=8&f<=50;dts=zeros(3,1);
P=zeros(N/2+1,4);
for w=1:4
 k=t>=win(w,1)&t<win(w,2);L=nnz(k);x=D(k,:).*(.5-.5*cos(2*pi*(0:L-1)'/(L-1)));F=abs(fft(x,N)).^2;P(:,w)=mean(F(1:N/2+1,:),2);
end
for w=2:4,y=log(P(band,w)./P(band,1));p=polyfit(f(band),y,1);dts(w-1)=-p(1)/(2*pi);end % power: exp(-2 pi f dt*)
[a0,~]=ray_absorption(t,offset,cm);ctr=mean(win,2);tstar=interp1(t,a0,ctr);dt1=tstar(2:4)-tstar(1);
c0=sum(dt1.^2)/sum(dt1.*dts);[~,i]=min(abs(log(grid)-log(max(c0,eps))));c=grid(i);
end
