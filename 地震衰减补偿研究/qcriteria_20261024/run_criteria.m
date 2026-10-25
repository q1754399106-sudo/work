function run_criteria()
% Non-sparse Q-selection criteria (inputs/PLAN_CRITERIA_ZH.md). Data-only selection;
% c_opt is read from the archived variant runs for scoring.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(2);
grid=exp(linspace(log(.6),log(2.2),25));traces=1:3:99;rows={};started=tic;
ref=[jsondecode(fileread(fullfile(root,'analysis','variants_F1_F2.json')));jsondecode(fileread(fullfile(root,'analysis','variants_F3_F4.json')));
 jsondecode(fileread(fullfile(root,'analysis','variants_F1_F2_F3_F4_M5_blocky_sparse.json')))];
for s={'F1','F2','F3','F4'}
 in=load(fullfile(root,'dev',[s{1} '.mat']));T=load(fullfile(root,'dev',[s{1} '_truth.mat']),'base');
 cm=in.config;cm.t=cm.t(:);offset=cm.solver_offset;sig=reshape(in.obs(1).sigma(traces),1,[]);
 for data={'M1_full','M5_blocky_sparse','M1_noisy'}
  if strcmp(data{1},'M1_noisy'),D=in.obs(1).D(:,traces);copt=ref(strcmp({ref.set},s{1})&strcmp({ref.variant},'M1_full')).c_opt;
  else,D=variant_data(T.base,cm,traces,data{1});copt=ref(strcmp({ref.set},s{1})&strcmp({ref.variant},data{1})).c_opt;end
  crit=zeros(4,numel(grid));
  for g=1:numel(grid)
   c=cm;c.Q=cm.Q*grid(g);d=prune_dictionary(c,offset);A=d.G*d.W;crit(:,g)=criteria(D,A,d,sig);
  end
  sel=zeros(1,5);for k=1:4,[~,i]=min(crit(k,:));sel(k)=grid(i);end
  sel(5)=spectral_factor(D,cm,offset,grid);
  names={'C0_QScan','C1_evidence','C2_GCV','C3_holdout','C4_spectral'};
  for k=1:5
   rows{end+1}=struct('set',s{1},'data',data{1},'criterion',names{k},'c_sel',sel(k),'c_opt',copt,'bias_log',log(sel(k)/copt)); %#ok<AGROW>
  end
  fprintf('%s %-17s c_opt %.3f | C0 %.3f C1 %.3f C2 %.3f C3 %.3f C4 %.3f | %.0fs\n',s{1},data{1},copt,sel,toc(started));
 end
end
uq_writejson(fullfile(root,'analysis','criteria.json'),[rows{:}]);
end

function v=criteria(D,A,d,sig)
% [C0 lasso objective; C1 Gaussian evidence; C2 GCV; C3 5-fold block hold-out]
[n,m]=size(D);lam=.5*sqrt(2*log(n))*sig;
X=refit_lasso(d.Gram,d.B'*D,lam,d.L);res=D-d.B*X;C0=sum((.5*sum(res.^2,1)+lam.*sum(abs(X),1))./sig.^2);
[U,S,~]=svd(A);s2=diag(S).^2;P=(U'*D).^2;s2grid=logspace(-8,2,81);C1=0;mu=zeros(1,m);
for j=1:m
 den=sig(j)^2+s2grid.*s2;ev=sum(log(den),1)+sum(P(:,j)./den,1);[e,i]=min(ev);C1=C1+.5*e;mu(j)=sig(j)^2/s2grid(i);
end
mugrid=logspace(-10,2,61)*max(s2);C2=0;
for j=1:m
 phi=s2./(s2+mugrid);r=sum(((1-phi).^2).*P(:,j),1);g=n*r./(n-sum(phi,1)).^2;C2=C2+min(g);
end
folds=round(linspace(0,n,6));C3=0;
for f=1:5
 te=folds(f)+1:folds(f+1);tr=setdiff(1:n,te);[Ut,St,Vt]=svd(A(tr,:),'econ');st=diag(St);B=Ut'*D(tr,:);
 for j=1:m,x=Vt*((st./(st.^2+mu(j))).*B(:,j));C3=C3+sum((A(te,:)*x-D(te,j)).^2)/sig(j)^2;end
end
v=[C0;C1;C2;C3];
end

function c=spectral_factor(D,cm,offset,grid)
% Spectral-ratio t* differences between time windows vs dictionary prediction.
win=[.30 .60;.60 .90;.90 1.20;1.20 1.60];t=cm.t;N=512;f=(0:N/2)'/(N*cm.dt);band=f>=8&f<=50;dts=zeros(3,1);
P=zeros(N/2+1,4);
for w=1:4
 k=t>=win(w,1)&t<win(w,2);L=nnz(k);x=D(k,:).*(.5-.5*cos(2*pi*(0:L-1)'/(L-1)));F=abs(fft(x,N)).^2;P(:,w)=mean(F(1:N/2+1,:),2);
end
for w=2:4,y=log(P(band,w)./P(band,1));p=polyfit(f(band),y,1);dts(w-1)=-p(1)/(2*pi);end % power: exp(-2 pi f dt*)
[a0,~]=ray_absorption(t,offset,cm);ctr=mean(win,2);tstar=interp1(t,a0,ctr);dt1=tstar(2:4)-tstar(1);
c0=sum(dt1.^2)/sum(dt1.*dts);[~,i]=min(abs(log(grid)-log(max(c0,eps))));c=grid(i);
end
