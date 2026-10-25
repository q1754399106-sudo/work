function test_coherence_math()
% Algebraic limits and a controlled coherent/isolated example.
n=200;dt=.004;f=30;t=(0:n-1)'*dt;z=pi*f*(t-t');W=(1-2*z.^2).*exp(-z.^2);
d=struct('W',W,'B',W./sqrt(sum(W.^2,1)),'scale',sqrt(sum(W.^2,1))','frequency',f,'dt',dt);
d.source_gram=(W./d.scale')'*(W./d.scale');
nx=5;R0=zeros(n,nx);supports=cell(1,nx);E=supports;V=supports;dicts=repmat({d},1,nx);
for j=1:nx
 s=80+2*j;R0(s,j)=.02*(1+.1*j);supports{j}=s;
end
R0(150,3)=.02;supports{3}=[supports{3};150];
R0(120,1)=-.03;supports{1}=sort([supports{1};120]);
R0(40,5)=.025;supports{5}=sort([40;supports{5}]);
for j=1:nx,ng=numel(unique(cumsum([1;diff(supports{j})>4])));E{j}=ones(ng,1);V{j}=ones(ng,1);end
live=true(1,nx);
base=struct('half_width',2,'dip_periods',.25,'lambda_min',1,'lambda_max',1);
same=coherence_soft_gather(R0,supports,dicts,E,V,live,base);
assert(all(cellfun(@(g)all(abs(g-.5)<1e-15),same.group_gain)),'lambda=1 must equal soft gain.');
opts=base;opts.lambda_min=.25;opts.lambda_max=4;
out=coherence_soft_gather(R0,supports,dicts,E,V,live,opts);
c3=out.coherence{3};assert(numel(c3)==2);
assert(c3(1)>.7,'Dipping event should be supported by aligned neighbours.');
assert(c3(2)<.05,'Isolated wavelet should not be supported.');
assert(out.group_gain{3}(1)>.5&&out.group_gain{3}(2)<.5);
assert(all(abs(out.reflectivity(:))<=abs(R0(:))+1e-15)&&all(sign(out.reflectivity(R0~=0))==sign(R0(R0~=0))));
flat=opts;flat.dip_periods=0;noalign=coherence_soft_gather(R0,supports,dicts,E,V,live,flat);
assert(noalign.coherence{3}(1)<c3(1),'Without lag alignment a dipping event loses support.');
dead=live;dead(4)=false;fb=coherence_soft_gather(R0,supports,dicts,E,V,dead,opts);
assert(all(fb.fallback{3})&&all(abs(fb.group_gain{3}-.5)<1e-15),'Near a bad trace the old rule applies.');
% A reflector that stops (fault/pinch-out) keeps one-sided support.
R1=R0;R1(:,4:5)=0;sup1=supports;sup1{4}=[];sup1{5}=[];E1=E;V1=V;E1{4}=[];V1{4}=[];E1{5}=[];V1{5}=[];
end1=coherence_soft_gather(R1,sup1,dicts,E1,V1,live,opts);
assert(end1.coherence{3}(1)>.7,'Terminating event should keep support from its continuous side.');
fprintf('Coherence gain limits, alignment, isolation, fallback and termination checks passed.\n');
end
