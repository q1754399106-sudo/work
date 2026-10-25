function c=cv_normalize_context(c)
active=sum(full(sum(c.H.^2,2))>0);s=norm(c.H,'fro')/sqrt(max(active,1));if s==0,s=1;end
c.H=c.H/s;c.HH=c.HH/s^2;
for j=1:c.nx,c.HHblocks{j}=c.HHblocks{j}/s^2;end
c.H_rms=s;
end
