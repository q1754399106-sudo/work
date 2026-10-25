function datasets=gen_layered_scene(cm,reflectivity,conditions)
% Four-layer synthesis as in q_dense_gathers run_dense/prepare_dense, new noise draws:
% AR(1)=0.45 temporal colouring plus 10% trace-mean (spatially common) component.
t=cm.t(:);z=pi*cm.fpeak*(t-t');W=(1-2*z.^2).*exp(-z.^2);target=W*reflectivity;
[a,b]=ray_absorption(t,450,cm);clean=attenuation_matrix(a+b,cm)*target;[nt,nx]=size(target);
datasets=struct('name',{},'snr_db',{},'seed',{},'D',{},'target',{},'sigma',{});
for k=1:numel(conditions)
 c=conditions(k);rng(c.seed,'twister');u=randn(nt,nx);
 for i=2:nt,u(i,:)=u(i,:)+.45*u(i-1,:);end
 u=.9*u+.1*mean(u,2);N=u*(norm(clean,'fro')/(10^(c.snr_db/20)*norm(u,'fro')));D=clean+N;
 sigma=zeros(1,nx);for j=1:nx,ni=estimate_trace_snr(D(:,j),t);sigma(j)=ni.sigma;end
 datasets(k)=struct('name',c.name,'snr_db',c.snr_db,'seed',c.seed,'D',D,'target',target,'sigma',sigma);
end
end
