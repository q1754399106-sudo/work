function prep_data(name)
% Datasets of model name in the frequency domain, cache/data_<name>.mat:
%   ref  elastic, no noise (scoring target only)        att  viscoacoustic, no noise
%   nXX  viscoacoustic + band-limited (3-55 Hz) noise, independent per trace, rms = rms(att in
%        1.2-1.8 s, |offset| <= 600 m) / 10^(XX/20); XX in P.noise_db (seed per level)
% Direct wave removed with the T0 records (same offset). Df: nrec x nf x nshot (single complex).
% Also stores the time-domain att/noisy records for route C.
P=rtm_setup();ns=numel(P.ixs);nr=numel(P.ixr);T0=load(fullfile(P.cache,'shots_T0','s000.mat'));
ATT=zeros(P.nt,nr,ns,'single');ELA=ATT;
for k=1:ns
 S=load(fullfile(P.cache,['shots_' name],sprintf('s%03d.mat',k)));
 o=abs(P.ixr-P.ixs(k));idx=o+3;                             % T0: source at interior column 3
 ATT(:,:,k)=S.att-T0.att(:,idx);ELA(:,:,k)=S.ela-T0.ela(:,idx);
end
off=(P.ixr(:)-P.ixs(:)')*P.dx;tw=P.t>=1.2&P.t<=1.8;near=abs(off)<=600;
a=ATT(tw,:,:);a=reshape(a,nnz(tw),[]);rms_s=sqrt(mean(a(:,near(:)).^2,'all'));
E=single(exp(-1i*2*pi*P.t*P.f));
D.ref=dft(ELA,E);D.att=dft(ATT,E);TD.att=ATT;
fq=(0:P.nt-1)'/(P.nt*P.dt);fq=min(fq,1/P.dt-fq);band=ones(P.nt,1);
band(fq<3|fq>55)=0;b=fq>=3&fq<6;band(b)=.5*(1-cos(pi*(fq(b)-3)/3));b=fq>50&fq<=55;band(b)=.5*(1+cos(pi*(fq(b)-50)/5));
for L=P.noise_db
 rng(20261117+round(100*(L+50)),'twister');n=randn(P.nt,nr*ns);n=real(ifft(fft(n).*band));
 n=single(n/sqrt(mean(n(:).^2))*rms_s/10^(L/20));X=ATT+reshape(n,P.nt,nr,ns);
 tag=sprintf('n%02d',L);if L<0,tag=sprintf('m%02d',-L);end
 D.(tag)=dft(X,E);TD.(tag)=X;
end
save(fullfile(P.cache,['data_' name '.mat']),'-struct','D','-v7.3');
save(fullfile(P.cache,['tdata_' name '.mat']),'-struct','TD','-v7.3');
fprintf('prep_data %s: signal rms (1.2-1.8 s, |o|<=600 m) %.3g\n',name,rms_s);
end

function Df=dft(X,E)
[nt,nr,ns]=size(X);Df=reshape(E.'*reshape(X,nt,nr*ns),[],nr,ns);Df=permute(Df,[2 1 3]);
end
