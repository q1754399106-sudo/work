function make_figures(name)
% Descriptive figures for stage 1: model, image panels, S-m curves (after score_stage1).
P=rtm_setup();I=load(fullfile(P.cache,['images_' name '.mat']));M=make_model(name,P);
R=jsondecode(fileread(fullfile(P.root,'analysis',['stage1_' name '.json'])));
out=fullfile(P.root,'figures');if ~isfolder(out),mkdir(out);end
z=P.z(P.in{1});x=P.x(P.in{2});set(0,'DefaultAxesFontSize',9);
% model
f=figure('Visible','off','Position',[0 0 1100 380]);
subplot(1,2,1);imagesc(x,z,M.v(P.in{1},P.in{2}));axis image;colorbar;title('velocity (m/s)');xlabel('x (m)');ylabel('z (m)');
subplot(1,2,2);imagesc(x,z,M.Q(P.in{1},P.in{2}));axis image;colorbar;title('Q');xlabel('x (m)');
exportgraphics(f,fullfile(out,['model_' name '.png']),'Resolution',150);close(f);
% images at m = 0.7 and 1 on m12
ds='m12';ms=I.meta.m_full{strcmp(I.meta.m_full(:,1),ds),2};
pan={'ref','R0',NaN;ds,'R0',NaN;ds,'IQ70',.7;ds,'W0',.7;ds,'UQW0',.7;ds,'IQ70',1;ds,'W0',1;ds,'UQW0',1};
f=figure('Visible','off','Position',[0 0 1400 700]);
for k=1:size(pan,1)
 [d,q,m]=pan{k,:};
 if strcmp(q,'R0'),img=I.(d).R0;tt=sprintf('%s R0',d);else,img=I.(d).(q)(:,ms==m);tt=sprintf('%s %s m=%.1f',d,q,m);end
 img=prep(img,P);c=prctile(abs(img(z>=820,:)),99,'all');
 subplot(2,4,k);imagesc(x,z,img,[-c c]);colormap(gray);axis image;title(tt,'Interpreter','none');
end
exportgraphics(f,fullfile(out,['images_' name '.png']),'Resolution',150);close(f);
% S-m curves
f=figure('Visible','off','Position',[0 0 1200 650]);ws={'deep','lens'};dss={'n00','m12'};k=0;
for a=1:2
 for b=1:2
  k=k+1;subplot(2,2,k);T=R.(dss{a}).(ws{b});hold on
  for q={'R0','PH','IQ30','IQ70','W0','UQW0','IQoracle'}
   y=T.(q{1});if isscalar(y),y=y*ones(size(T.m));end
   lw=1;if strcmp(q{1},'UQW0'),lw=2.2;end
   semilogx(T.m,y,'-o','LineWidth',lw,'MarkerSize',3,'DisplayName',q{1});
  end
  set(gca,'XScale','log');xticks(T.m);grid on;title(sprintf('%s, W_%s',dss{a},ws{b}),'Interpreter','none');xlabel('Q multiplier m');ylabel('S')
  if k==1,legend('Location','southeast');end
 end
end
exportgraphics(f,fullfile(out,['S_m_' name '.png']),'Resolution',150);close(f);
end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=10;h=exp(-.5*((-30:30)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),30,1);y;repmat(y(end,:),30,1)],'valid');
y(P.z(P.in{1})<140,:)=0;
end
