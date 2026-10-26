function prior=uq_build_prior(D,models,c,domain,preset)
assert(any(strcmp(domain,{'Offset','Section'})));
[n,nx]=size(D); sigma=zeros(1,nx); lambda=sigma; pilot=zeros(n,nx);
for j=1:nx
 ni=estimate_trace_snr(D(:,j),c.t); sigma(j)=ni.sigma;
 [~,lambda(j)]=discrepancy_apply(models{j}.direct,D(:,j),c.lambda_grid,sqrt(n)*sigma(j));
 F=models{j}.direct; u=F.U'*D(:,j); pilot(:,j)=F.V*(F.z./(F.z.^2+preset.pilot_multiplier*lambda(j)).*u);
end
H=sparse(n*nx,n*nx);
if preset.gamma>0
 [H,info]=learn_structure_operator(pilot,c,structure_options());
 if strcmp(domain,'Offset'), R=avo_projector_v2(models); H=H*R; end
end
prior=struct('H',H,'sigma',sigma,'lambda',lambda,'pilot',pilot);
end
