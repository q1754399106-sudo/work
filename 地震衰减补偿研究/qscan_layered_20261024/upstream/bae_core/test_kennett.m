function info=test_kennett()
% (1) Elastic single interface: Ricker peak equals r at 2h/v.
% (2) Attenuated single interface vs the dictionary forward model (same Q):
%     the mismatch (Kjartansson vs Kolsky-Futterman) is reported, bounded loosely.
cm=struct('dt',.004,'nt',451,'nfft',4096,'fref',100,'fpeak',32);cm.t=(0:cm.nt-1)'*cm.dt;
v=[2000 2400];rho=[2.1 2.2];h=600;Q=[50 50];
[tr,el]=kennett_trace(v,rho,h,Q,cm,32);r=(v(2)*rho(2)-v(1)*rho(1))/(v(2)*rho(2)+v(1)*rho(1));
[pk,i]=max(el);t0=2*h/v(1);
assert(abs(pk-r)<.01*abs(r)&&abs(cm.t(i)-t0)<=cm.dt,'Elastic single-interface check failed.');
c=struct('dt',cm.dt,'nt',cm.nt,'nfft',cm.nfft,'fref',cm.fref,'vp',[2000 2000 2000 2000],'h',[300 300 300],'Q',[50 50 50 50],'fpeak',32,'t',cm.t);
d=prune_dictionary(c,0);k=round(t0/cm.dt)+1;e=zeros(cm.nt,1);e(k)=r;model=d.G*d.W*e;
mismatch=norm(model-tr)/norm(tr);assert(mismatch<.2,'Attenuated forward mismatch too large: %.3f',mismatch);
info=struct('elastic_peak',pk,'r',r,'attenuated_vs_dictionary_relative_difference',mismatch);
fprintf('Kennett checks passed: elastic peak %.4f vs r %.4f; attenuated vs dictionary model %.3f\n',pk,r,mismatch);
end
