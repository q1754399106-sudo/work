function R=score_block(names)
% Step 1 of PLAN_BLOCKFLIP_ZH.md: accuracy of the block-polarity 1/snr against the truth in W_deep.
% e = median over f in [10, 50] Hz of the median over W_deep pixels of log2(est / true).
% Verdicts A1, A2 and the choice of B*; writes analysis/blockflip_A.json.
if nargin<1,names={'D1','D2','D3','D4'};end
P=rtm_setup();z=P.z(P.in{1});x=P.x(P.in{2});wd=(z>=820&z<=1550)'&(x>=300&x<=2100);wd=wd(:);
sets={'m30','m40','g30','l30'};fb=P.f>=10&P.f<=50;fr=[10 25 40];
R=struct('BS',[],'models',{names});E=struct();
for n=names
 I=load(fullfile(P.cache,['invblk_' n{1} '.mat']));BS=I.BS;R.BS=BS;
 for s=sets
  tv=I.([s{1} '_true'])(wd,:);
  for B=BS
   lr=log2(double(I.(sprintf('%s_B%d',s{1},B))(wd,:))./double(tv));mf=median(lr,1);
   T=struct('e',median(mf(fb)),'e_f',mf,'at_10_25_40',mf(ismember(P.f,fr)), ...
            'median_true_10_25_40',median(tv(:,ismember(P.f,fr)),1));
   E.(n{1}).(s{1}).(sprintf('B%d',B))=T;
  end
 end
end
R.e=E;R.f=P.f;
main={'m30','m40','l30'};cand=[8 32 300];worst=zeros(size(cand));
for c=1:numel(cand)
 w=0;for n=names,for s=main,w=max(w,abs(E.(n{1}).(s{1}).(sprintf('B%d',cand(c))).e));end,end;worst(c)=w;
end
[~,ib]=min(worst);Bs=cand(ib);R.worst_abs_e=struct('B',cand,'value',worst);R.Bstar=Bs;
A1=true;A2=true;
for n=names
 A1=A1&&abs(E.(n{1}).l30.(sprintf('B%d',Bs)).e)<=1;
 for s={'m30','m40'},A2=A2&&abs(E.(n{1}).(s{1}).(sprintf('B%d',Bs)).e)<=abs(E.(n{1}).(s{1}).B1.e)+.5;end
end
R.A1=A1;R.A2=A2;R.pass=A1&&A2;
fid=fopen(fullfile(P.root,'analysis','blockflip_A.json'),'w');fprintf(fid,'%s',jsonencode(R,'PrettyPrint',true));fclose(fid);
fprintf('B* = %d (worst |e| %s for B = %s); A1 %d A2 %d\n',Bs,mat2str(worst,3),mat2str(cand),A1,A2);
for n=names
 for s=sets
  fprintf('%s %s e:',n{1},s{1});for B=R.BS,fprintf(' B%d %+.2f',B,E.(n{1}).(s{1}).(sprintf('B%d',B)).e);end;fprintf('\n');
 end
end
end
