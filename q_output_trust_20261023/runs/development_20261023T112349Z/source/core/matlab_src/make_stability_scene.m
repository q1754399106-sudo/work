function scene=make_stability_scene(name,c)
% Synthetic geometry only. No true scene geometry is passed to the inverse.
nx=11; offsets=0:150:1500; lateral=(0:nx-1)*50;
if endsWith(name,'Section'), offsets=750*ones(1,nx); end
n=c.nt; truth=zeros(n,nx); events=cell(1,nx); coeff=events; templates=events; ang=events;
for j=1:nx
 [~,~,a]=make_reference(offsets(j),c); a=a(:); times=c.event_t(:);
 if strcmp(name,'FaultSection')
  times=times+.008*sin((j-1)*pi/10); if j>=7, times(2:3)=times(2:3)+.04; end
 elseif strcmp(name,'PinchoutSection')
  times=[c.event_t(1);c.event_t(2);c.event_t(2)+.030*max(0,(9-j)/8)]; a=[a(1);a(2);-a(2)];
 elseif strcmp(name,'WeakSection')
  times(2)=times(2)+.012*sin((j-1)*pi/10); a(2)=.08*a(2);
 end
 W=zeros(n,3);
 for k=1:3, z=pi*c.fpeak*(c.t-times(k)); W(:,k)=(1-2*z.^2).*exp(-z.^2); end
 truth(:,j)=W*a; events{j}=times; coeff{j}=a; templates{j}=W;
 [~,~,~,p]=ray_absorption(c.event_t,offsets(j),c); ang{j}=(p(:).*c.vp(1:3)').^2;
end
scene=struct('name',name,'truth',truth,'offsets',offsets,'lateral',lateral,'events',{events},'coeff',{coeff},'templates',{templates},'sin2',{ang});
end
