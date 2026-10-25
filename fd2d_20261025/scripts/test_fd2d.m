function test_fd2d()
% Verification of the 2-D forward and its processing, on test geologies only (none of H1-H6):
%  T1 staggered 8th-order derivatives (conv2 alignment) against the explicit stencil;
%  T2 homogeneous viscoacoustic medium: attenuated/elastic spectral ratio at 400 m against the
%     exact 2-D solution (Hankel functions with the GSLS modulus), and the absolute elastic
%     response against the 2-D Green's function;
%  T3 laterally homogeneous fine-layered column (carbonate, gas sand): processed zero-offset
%     trace against kennett_trace (exact 1-D plane-wave response, Kjartansson Q);
%  T4 absorbing boundaries: one flat reflector, +-500 m window against +-1000 m window;
%  T5 GSLS: achieved Q(f) on 5-100 Hz for Q = 18, 40, 150.
% Writes analysis/verification.json and cache/verify/T3.mat.
P=fd_setup();maxNumCompThreads(P.threads);fd=P.fd;out=struct();
%% T1
c=[1225/1024 -245/3072 49/5120 -5/7168];w=[-fliplr(c) c];A=randn(40,50);
D=conv2(A,fliplr(w),'same');Dm=conv2(A,fliplr([w 0]),'same');E=zeros(40,50);Em=E;
for j=5:46,for k=1:4,E(:,j)=E(:,j)+c(k)*(A(:,j+k)-A(:,j+1-k));Em(:,j)=Em(:,j)+c(k)*(A(:,j+k-1)-A(:,j-k));end,end
out.T1=struct('forward_max_err',max(abs(D(:,5:46)-E(:,5:46)),[],'all'),'backward_max_err',max(abs(Dm(:,5:46)-Em(:,5:46)),[],'all'));
fprintf('T1 derivative alignment: %.2e %.2e\n',out.T1.forward_max_err,out.T1.backward_max_err);
%% T2
v0=2000;r0=.31*v0^.25;h=round(1000/fd.dx)+fd.npml;n=2*h+1;nt=round(.8/fd.dt);
M=struct('v',v0*ones(n),'rho',r0*ones(n),'Q',40*ones(n));
ir=round(400/fd.dx);opt=struct('dx',fd.dx,'dt',fd.dt,'nt',nt,'fref',P.cm.fref,'npml',fd.npml,'src',[h+1 h+1],'rec',[h+1 h+1+ir;h+1+ir h+1]);
tq=(1:nt)'*fd.dt-fd.t0;a=(pi*P.cm.fpeak*tq).^2;q=(1-2*a).*exp(-a);
Ra=fd2d_visco(M,q,opt);Me=M;Me.Q(:)=Inf;Re=fd2d_visco(Me,q,opt);
N=4096;f=(0:N/2)'/(N*fd.dt);wv=2*pi*f;band=f>=8&f<=90;r=ir*fd.dx;
tw=ones(nt,1);k0=round(.7/fd.dt);tw(k0:end)=.5*(1+cos(pi*(0:nt-k0)'/(nt-k0)));   % taper before boundary returns
Fa=fft(Ra(:,1).*tw,N);Fe=fft(Re(:,1).*tw,N);Fa=Fa(1:N/2+1);Fe=Fe(1:N/2+1);
G=gsls_design(40,P.cm.fref);m=1+G.tau*sum(1i*wv.*G.ts./(1+1i*wv.*G.ts),2);MR=r0*v0^2*G.vfac;
ka=wv.*sqrt(r0./(MR*m));ka=real(ka)-1i*abs(imag(ka));ke=wv/v0;
Rth=(r0*v0^2./(MR*m)).*besselh(0,2,ka*r)./besselh(0,2,ke*r);Rfd=Fa./Fe;e=Rfd(band)./Rth(band)-1;
Qs=fft(q,N)*fd.dt;Qs=Qs(1:N/2+1);Gabs=(wv*r0/(4*r0*v0^2)).*Qs.*besselh(0,2,ke*r);ea=(Fe(band)*fd.dt)./Gabs(band)-1;
vert=fft(Ra(:,2).*tw,N);vert=vert(1:N/2+1);
out.T2=struct('ratio_max_abs_err',max(abs(e)),'ratio_rms_err',rms(abs(e)),'ratio_max_phase_deg',max(abs(angle(1+e)))*180/pi, ...
 'abs_elastic_median_amp_err',median(abs(abs(1+ea)-1)),'abs_elastic_median_phase_deg',median(angle(1+ea))*180/pi, ...
 'isotropy_max_rel',max(abs(vert(band)-Fa(band)))/max(abs(Fa(band))));
fprintf('T2 homogeneous Q=40 at 400 m: ratio err max %.3f rms %.3f, phase max %.2f deg | absolute elastic amp %.3f phase %.2f deg | isotropy %.3f\n', ...
 out.T2.ratio_max_abs_err,out.T2.ratio_rms_err,out.T2.ratio_max_phase_deg,out.T2.abs_elastic_median_amp_err,out.T2.abs_elastic_median_phase_deg,out.T2.isotropy_max_rel);
%% T3
sp=struct('name','T','seed',2026102799,'macro_vp',[1800 2300 2800 3300],'macro_h',[230 300 350],'macro_Q',[40 60 85 110], ...
 'carbonate',true,'fault_trace',1e6,'throw',0);
Gm=fd_geology(sp,P);[~,ix50]=min(abs(Gm.x-Gm.trace_x(50)));
for f_={'v','rho','Q'},Gm.(f_{1})=repmat(Gm.(f_{1})(:,ix50),1,numel(Gm.x));end
[a3,e3]=fd_shot(Gm,ix50,P);Hm=homog(Gm);[ah,eh]=fd_shot(Hm,ix50,P,round(fd.r0/fd.dx));
cfg=Gm.config;fa=process_zo(a3,ah(:,1),eh(:,2),Gm.cap_v,cfg,P);fe=process_zo(e3,eh(:,1),eh(:,2),Gm.cap_v,cfg,P);
c5=Gm.col50;[ka3,ke3]=kennett_trace(c5.v,c5.rho,c5.h,c5.Q,P.cm,P.cm.fpeak);kg3=kennett_gsls(c5.v,c5.rho,c5.h,c5.Q,P.cm,P.cm.fpeak);
win=P.cm.t>=.25;late=P.cm.t>=1;rel=@(x,y) norm(x(win)-y(win))/norm(y(win));cc=@(x,y) corr(x(win),y(win));
out.T3=struct('rel_err_attenuated',rel(fa,ka3),'rel_err_elastic',rel(fe,ke3),'corr_attenuated',cc(fa,ka3),'corr_elastic',cc(fe,ke3), ...
 'energy_ratio_fd',norm(fa(win))/norm(fe(win)),'energy_ratio_kennett',norm(ka3(win))/norm(ke3(win)), ...
 'late_rel_err_attenuated',rel2(fa,ka3,late),'late_rel_err_elastic',rel2(fe,ke3,late), ...
 'rel_err_attenuated_vs_kennett_gsls',rel(fa,kg3),'late_rel_err_attenuated_vs_kennett_gsls',rel2(fa,kg3,late), ...
 'late_rel_err_kennett_gsls_vs_kjartansson',rel2(kg3,ka3,late));
fprintf(['T3 layered vs Kennett: rel err att %.3f ela %.3f | corr att %.3f ela %.3f | att/ela energy FD %.3f Kennett %.3f | t>=1 s rel err att %.3f ela %.3f\n' ...
 '   vs 1-D GSLS response: att %.3f, t>=1 s %.3f | 1-D GSLS vs Kjartansson, t>=1 s %.3f\n'], ...
 out.T3.rel_err_attenuated,out.T3.rel_err_elastic,out.T3.corr_attenuated,out.T3.corr_elastic,out.T3.energy_ratio_fd,out.T3.energy_ratio_kennett, ...
 out.T3.late_rel_err_attenuated,out.T3.late_rel_err_elastic,out.T3.rel_err_attenuated_vs_kennett_gsls,out.T3.late_rel_err_attenuated_vs_kennett_gsls, ...
 out.T3.late_rel_err_kennett_gsls_vs_kjartansson);
vdir=fullfile(P.root,'cache','verify');if ~isfolder(vdir),mkdir(vdir);end
save(fullfile(vdir,'T3.mat'),'fa','fe','ka3','ke3','kg3','cfg','-v7');
%% T4
G4=Gm;G4.v(:)=2000;G4.Q(:)=Inf;G4.v(G4.z>=800,:)=2600;G4.rho=.31*G4.v.^.25;
P2=P;P2.fd.half_window=1000;H4=G4;H4.v(:)=2000;H4.rho(:)=.31*2000^.25;
[~,e4]=fd_shot(G4,ix50,P);[~,h4]=fd_shot(H4,ix50,P,round(fd.r0/fd.dx));
G4w=widen(G4,P2);H4w=widen(H4,P2);[~,ix]=min(abs(G4w.x-Gm.trace_x(50)));
[~,e4w]=fd_shot(G4w,ix,P2);[~,h4w]=fd_shot(H4w,ix,P2,round(fd.r0/fd.dx));
cfg4=cfg;cfg4.vp=[2000 2600 2600 2600];cfg4.h=[800 500 500];   % spreading from this model's own velocities
t4=process_zo(e4,h4(:,1),h4(:,2),2000,cfg4,P);t4w=process_zo(e4w,h4w(:,1),h4w(:,2),2000,cfg4,P2);
Rn=(2600^1.25-2000^1.25)/(2600^1.25+2000^1.25);
Nf=1024;fq=(0:Nf/2)'/(Nf*P.cm.dt);S=fft(t4w,Nf)*P.cm.dt;S=S(1:Nf/2+1);src=(2/sqrt(pi))*(fq.^2/P.cm.fpeak^3).*exp(-(fq/P.cm.fpeak).^2);
b=fq>=15&fq<=50;amp=median(abs(S(b))./(Rn*src(b)));
out.T4=struct('window_rel_diff',norm(t4(win)-t4w(win))/norm(t4w(win)),'amplitude_over_R_15_50Hz',amp);
fprintf('T4 boundary: +-500 vs +-1000 m rel diff %.4f; reflection amplitude / R (15-50 Hz spectral ratio) %.3f\n',out.T4.window_rel_diff,amp);
%% T5
G5=gsls_design([18;40;150],P.cm.fref);k=G5.f>=5&G5.f<=100;dev=max(abs(G5.q(:,k)./[18;40;150]-1),[],2);
out.T5=struct('Q',[18 40 150],'max_rel_dev_5_100Hz',dev');
fprintf('T5 GSLS Q(f) max deviation on 5-100 Hz: %s\n',mat2str(dev',3));
A_=fullfile(P.root,'analysis');if ~isfolder(A_),mkdir(A_);end
fid=fopen(fullfile(A_,'verification.json'),'w');fprintf(fid,'%s',jsonencode(out,'PrettyPrint',true));fclose(fid);
end

function H=homog(Gm)
H=Gm;H.v(:)=Gm.cap_v;H.rho(:)=Gm.cap_rho;H.Q(:)=Gm.cap_Q;
end

function W=widen(Gm,P2)
% Extend the model sideways by repeating the edge columns (+-1000 m windows).
add=round((P2.fd.half_window-500)/P2.fd.dx);W=Gm;
for f_={'v','rho','Q'},A=Gm.(f_{1});W.(f_{1})=[repmat(A(:,1),1,add) A repmat(A(:,end),1,add)];end
W.x=[Gm.x(1)-(add:-1:1)*P2.fd.dx Gm.x Gm.x(end)+(1:add)*P2.fd.dx];
end

function e=rel2(x,y,k)
e=norm(x(k)-y(k))/norm(y(k));
end

function r=corr(x,y)
c=corrcoef(x,y);r=c(1,2);
end

function tr=kennett_gsls(v,rho,h,Q,cm,fpeak)
% kennett_trace (bae_q_20261024/core) with the Kjartansson velocity replaced by the GSLS
% modulus of fd2d_visco: c(w) = v*sqrt(vfac*m(w)). Exact 1-D plane-wave response of the
% rheology the finite-difference code implements; isolates numerical and 2-D effects.
N=cm.nfft;f=(0:N/2)'/(N*cm.dt);w=2*pi*f;
src=(2/sqrt(pi))*(f.^2/fpeak^3).*exp(-(f/fpeak).^2);
[uq,~,iq]=unique(Q(:));G=gsls_design(uq,cm.fref);L=numel(v);C=zeros(numel(f),L);
for k=1:L
 g=iq(k);m=1+G.tau(g)*sum(1i*w.*G.ts./(1+1i*w.*G.ts),2);C(:,k)=v(k)*sqrt(G.vfac(g)*m);
end
Rd=zeros(numel(f),1);
for k=L:-1:2
 zk=rho(k)*C(:,k);zu=rho(k-1)*C(:,k-1);r=(zk-zu)./(zk+zu);
 if k<L,ph=exp(-2i*w.*h(k)./C(:,k));Rd=(r+Rd.*ph)./(1+r.*Rd.*ph);else,Rd=r;end
end
Rd=Rd.*exp(-2i*w.*h(1)./C(:,1));
X=src.*Rd;X(1)=0;X(end)=real(X(end));X=[X;conj(X(end-1:-1:2))];
tr=real(ifft(X))/cm.dt;tr=tr(1:cm.nt);
end
