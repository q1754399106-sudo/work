function test_bae_math()
% (1) zero model error reproduces FullCalRefit; (2) whitening identity.
cm=struct('dt',.004,'nt',201,'nfft',1024,'fref',100,'vp',[1800 2300 2800 3300],'h',[230 300 350],'Q',[40 60 85 110],'fpeak',32);
cm.t=(0:cm.nt-1)'*cm.dt;d=prune_dictionary(cm,450);rng(7,'twister');r=zeros(cm.nt,3);r([40 90 150],:)=[.1 -.05 .03]'*[1 1 1];
D=d.G*d.W*r+.002*randn(cm.nt,3);sig=[.002 .002 .002];initial=reflectivity_refit_prepared(D,repmat({d},1,3),sig);
same=bae_q_refit(D,d,repmat({d.G},1,4),sig,initial,true);
assert(norm(same.gather-initial.gather,'fro')<=1e-8*norm(initial.gather,'fro'),'Zero model error must reproduce FullCalRefit.');
E=randn(cm.nt,5);[U,Z]=svd(E,'econ');ev=diag(Z).^2;s=.3;shrink=1-1./sqrt(1+ev/s^2);
Wt=(eye(cm.nt)-U*diag(shrink)*U')/s;C=s^2*eye(cm.nt)+E*E';
assert(norm(Wt*C*Wt'-eye(cm.nt),'fro')<1e-9*cm.nt,'Whitening identity failed.');
fprintf('BAE math checks passed.\n');
end
