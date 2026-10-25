function [s,r,coeff] = make_reference(x,c)
% Ideal NMO-aligned gather. Shuey two-term AVO, not exact Zoeppritz.
% No geometric spreading, transmission, multiples, or NMO stretch.
[~,~,~,p]=ray_absorption(c.event_t,x,c);
r=zeros(c.nt,1); coeff=zeros(3,1); s=r;
for j=1:3
    vp=mean(c.vp(j:j+1)); vs=mean(c.vs(j:j+1)); rho=mean(c.rho(j:j+1));
    dvp=diff(c.vp(j:j+1))/vp; dvs=diff(c.vs(j:j+1))/vs;
    drho=diff(c.rho(j:j+1))/rho;
    A=0.5*(dvp+drho); B=0.5*dvp-2*(vs/vp)^2*(drho+2*dvs);
    coeff(j)=A+B*(p(j)*c.vp(j))^2;
    z=pi*c.fpeak*(c.t-c.event_t(j));
    s=s+coeff(j)*(1-2*z.^2).*exp(-z.^2);
    r(round(c.event_t(j)/c.dt)+1)=coeff(j);
end
end
