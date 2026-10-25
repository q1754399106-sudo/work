function out=reflectivity_soft_prepared(D,dicts,sigma,initial)
% Continuous source-domain noise-risk proxy on fixed FullCalRefit groups.
% No truth, model name, spatial target or user event locations are inputs.
% Noise propagation is conditional on the selected support and white noise.
[n,nx]=size(D);sigma=reshape(sigma,1,[]);S=zeros(n,nx);R=S;
gains=cell(1,nx);energies=gains;variances=gains;normals=zeros(1,nx);
for j=1:nx
 d=dicts{j};support=initial.support{j}(:);
 if isempty(support)
  gains{j}=[];energies{j}=[];variances{j}=[];continue;
 end
 C=d.B(:,support);[~,Z,V]=svd(C,'econ');z=diag(Z);
 top=max(z.^2);bottom=min(z.^2);
 reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));
 % K K' without constructing the nt-by-nt prediction matrix.
 noise_shape=(V.*reshape(z./(z.^2+reg),1,[]));
 noise_shape=noise_shape*noise_shape';
 a=initial.reflectivity(support,j).*d.scale(support);
 gap=max(1,round(.5/(d.frequency*d.dt)));
 group=cumsum([1;diff(support)>gap]);ng=max(group);
 gain=zeros(ng,1);energy=gain;variance=gain;weighted=a;
 for k=1:ng
  ii=find(group==k);Wg=d.source_gram(support(ii),support(ii));ag=a(ii);
  energy(k)=max(real(ag'*Wg*ag),0);
  variance(k)=max(real(sum(sum(Wg.*noise_shape(ii,ii)')))*sigma(j)^2,0);
  if variance(k)==0,gain(k)=1;
  else,gain(k)=energy(k)/(energy(k)+variance(k));end
  weighted(ii)=gain(k)*ag;
 end
 assert(all(isfinite(gain))&&all(gain>=0)&all(gain<=1));
 r=zeros(n,1);r(support)=weighted./d.scale(support);R(:,j)=r;S(:,j)=d.W*r;
 % Verify the *initial refit* normal equations. Shrinkage is a later step.
 normals(j)=norm(C'*(C*a-D(:,j))+reg*a)/max(norm(C'*D(:,j)),realmin);
 gains{j}=gain;energies{j}=energy;variances{j}=variance;
end
assert(all(isfinite(S(:)))&&max(normals)<1e-8);
out=initial;out.gather=S;out.reflectivity=R;out.group_gain=gains;
out.group_energy=energies;out.group_noise_energy=variances;
out.initial_normal_residual=max(normals);
out.selected_method='SoftGroupRefit';
out.rule='gain=source group energy/(source group energy+conditional propagated white noise energy)';
out.uncertainty_assumption='Fixed selected support; white noise proxy, no Q-model bias or selection uncertainty.';
end
