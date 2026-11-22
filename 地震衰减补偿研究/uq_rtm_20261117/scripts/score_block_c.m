function E=score_block_c(names)
% Descriptive accuracy of the 1/snr estimates (B = 1, 300) on the confirmation models (PLAN_BLOCKFLIP_CONFIRM_ZH.md),
% same metric as score_block.m: e = median over 10-50 Hz of the median over W_deep pixels of log2(est / true).
% Writes analysis/blockflip_C_accuracy.json.
if nargin<1,names={'C5','C6','C7','C8'};end
P=rtm_setup();z=P.z(P.in{1});x=P.x(P.in{2});wd=(z>=820&z<=1550)'&(x>=300&x<=2100);wd=wd(:);fb=P.f>=10&P.f<=50;E=struct();
for n=names
 I=load(fullfile(P.cache,['invblk_' n{1} '.mat']));
 for s={'m30','m40','l30'}
  tv=I.([s{1} '_true'])(wd,:);
  for B=I.BS
   mf=median(log2(double(I.(sprintf('%s_B%d',s{1},B))(wd,:))./double(tv)),1);E.(n{1}).(s{1}).(sprintf('B%d',B))=median(mf(fb));
  end
  fprintf('%s %s e: B1 %+.2f B300 %+.2f\n',n{1},s{1},E.(n{1}).(s{1}).B1,E.(n{1}).(s{1}).B300);
 end
end
fid=fopen(fullfile(P.root,'analysis','blockflip_C_accuracy.json'),'w');fprintf(fid,'%s',jsonencode(E,'PrettyPrint',true));fclose(fid);
end
