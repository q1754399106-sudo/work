function [D_work,b_coh,ortho_ratio,tail_mask,extra_diag]=p1_noise_uncertainty(D,pilot,models,sigma,lam,domain,R_avo,p1)
% Separate coherent interference handling in Offset domain and local noise-floor protection in Section domain.
[n,nx]=size(D); b_coh=zeros(n,nx); ortho_ratio=0;
tail_mask=false(n,nx); extra_diag=zeros(n,nx);
if strcmp(domain,'Offset')
 pilot_ortho=reshape(R_avo*pilot(:),n,nx);
 D_ortho=zeros(n,nx);
 for j=1:nx, D_ortho(:,j)=models{j}.G*pilot_ortho(:,j); end
 expected_ortho=((nx-2)/nx)*n*sum(sigma.^2);
 ortho_ratio=sum(D_ortho(:).^2)/max(expected_ortho,realmin);
 if ortho_ratio>=p1.coherent_ortho_ratio_gate
  lags=reshape(p1.coherent_prediction_lags,1,[]);
  t_idx=(1:n)';
  for j=1:nx
   neighbors=[j-1,j+1]; neighbors=neighbors(neighbors>=1 & neighbors<=nx);
   B=zeros(n,numel(neighbors)*numel(lags)); col=0;
   for jj=neighbors
    for lag=lags
     col=col+1; tt=min(n,max(1,t_idx+lag)); B(:,col)=D_ortho(tt,jj);
    end
   end
   pcol=size(B,2); ridge=p1.coherent_prediction_ridge*max(sum(B(:).^2)/pcol,realmin);
   coef=(B'*B+ridge*eye(pcol))\(B'*D_ortho(:,j));
   pred=B*coef;
   r2=1-sum((D_ortho(:,j)-pred).^2)/max(sum(D_ortho(:,j).^2),realmin);
   gate=min(1,max(0,(r2-p1.coherent_r2_gate)/p1.coherent_r2_span));
   b_coh(:,j)=gate*pred;
  end
 end
else
 k=ones(19,1)/19; loc_pow=conv2(pilot.^2,k,'same');
 lat_med=median(loc_pow,2);
 sigma_p2=zeros(1,nx);
 for j=1:nx
  z=models{j}.direct.z;
  sigma_p2(j)=(sigma(j)^2)*mean((z./(z.^2+4*lam(j))).^2);
 end
 tail_mask=(loc_pow<p1.section_tail_power_factor*sigma_p2) & (loc_pow<p1.section_tail_lateral_ratio*lat_med);
 for j=1:nx
  extra_diag(tail_mask(:,j),j)=p1.section_tail_diag_lambda_mult*lam(j);
 end
end
D_work=D-b_coh;
end
