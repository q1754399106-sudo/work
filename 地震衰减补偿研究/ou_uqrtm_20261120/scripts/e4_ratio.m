function e4_ratio()
% Step 0: two-way t* of the E4 / E2 compensation models divided by the t* of the correct model,
% over the deep window W_deep (all shots, specular receiver). Writes analysis/step0_tstar_ratio.json.
P=ou_setup();maxNumCompThreads(P.threads);P0=P;P0.cache=P.src_cache;
zi=P.z(P.in{1});xi=P.x(P.in{2});wd=(zi>=820&zi<=1550)'&(xi>=300&xi<=2100);J=struct();
for n={'D1','D2','D3','D4','C1','C2','C3','C4'}
 F0=load(fullfile(P.src_cache,['fields_' n{1} '.mat']));
 for key={'E4_070','E4_140','E2'}
  F1=load(fullfile(P.src_cache,['fields_' n{1} '_' key{1} '.mat']));r=[];
  for s=1:5:numel(P.ixs)
   [~,t0]=two_way(F0,s,P0);[~,t1]=two_way(F1,s,P0);q=t1(P.in{1},P.in{2})./t0(P.in{1},P.in{2});r=[r;q(wd)]; %#ok<AGROW>
  end
  J.(n{1}).(['x' key{1}])=prctile(r,[5 50 95]);
 end
 fprintf('%s E4_070 %s | E4_140 %s | E2 %s\n',n{1},mat2str(J.(n{1}).xE4_070,3),mat2str(J.(n{1}).xE4_140,3),mat2str(J.(n{1}).xE2,3));
end
fid=fopen(fullfile(P.root,'analysis','step0_tstar_ratio.json'),'w');fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true));fclose(fid);
end
