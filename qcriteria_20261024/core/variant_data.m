function [D0,Tg]=variant_data(base,cm,traces,variant)
% Noise-free observation and elastic target for 'M1_full' or 'M5_blocky_sparse',
% identical to run_variants.m.
nx=99;D0=zeros(cm.nt,numel(traces));Tg=D0;
for q=1:numel(traces)
 [v,rho,h,Q]=lateral(base,traces(q),nx);
 switch variant
  case 'M1_full',mode='full';
  case 'M5_blocky_sparse',mode='noloss';[v,rho,h,Q]=blocky(base,traces(q),nx);
 end
 [D0(:,q),Tg(:,q)]=kennett_variant(v,rho,h,Q,cm,cm.fpeak,mode);
end
end

function [v,rho,h,Q]=lateral(base,j,nx)
x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];hs=base.h.*stretch(base.macro);
hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);keep=hs>0;
v=[base.cap_v base.v(keep)];Q=[base.cap_Q base.Q(keep)];h=[base.cap_h hs(keep)];
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;h=h(1:end);
end

function [v,rho,h,Q]=blocky(base,j,nx)
[vf,~,hf,Qf]=lateral(base,j,nx);x=(j-1)/(nx-1);stretch=[1+.06*x 1+.05*sin(2*pi*x) 1-.05*x 1];
hs=base.h.*stretch(base.macro);hs(base.gas)=base.h(base.gas)*min(max((90-j)/50,0),1);m=base.macro(hs>0);
hh=hf(2:end);vv=vf(2:numel(hf));qq=Qf(2:numel(hf));v=vf(1);h=hf(1);Q=Qf(1);
for u=unique(m)
 k=m==u;tt=sum(hh(k)./vv(k));v(end+1)=sum(hh(k))/tt;h(end+1)=sum(hh(k));Q(end+1)=tt/sum(hh(k)./vv(k)./qq(k)); %#ok<AGROW>
end
v(end+1)=v(end)*1.05;Q(end+1)=Q(end);rho=.31*v.^.25;
end
