function gates()
% Registered gates (PLAN_ZH.md 1b H1-H4) for QScanRefit, plus descriptive means for all frozen methods.
root=fileparts(mfilename('fullpath'));out=struct();
for sp={'val','ind'}
 S=readtable(fullfile(root,'analysis',[sp{1} '_ratios_all.csv']),'TextType','char');methods=unique(S.method);res=struct();
 for q=1:numel(methods)
  r=S(strcmp(S.method,methods{q}),:);w=r.weak_q0_ratio(isfinite(r.weak_q0_ratio));
  m=struct('combos',height(r),'mean_wc_ratio',mean(r.wc_ratio),'mean_q0_ratio',mean(r.q0_ratio),'wc_improved',nnz(r.wc_ratio<1), ...
   'mean_weak_q0_ratio',mean(w),'worst_q0_ratio',max(r.q0_ratio));
  if strcmp(methods{q},'QScanRefit')
   m.H1_wc=m.mean_wc_ratio<=.95;m.H2_q0=m.mean_q0_ratio<=1.02;m.H3_majority=m.wc_improved>height(r)/2;
   m.H4_weak=isempty(w)||m.mean_weak_q0_ratio<=1.10;m.all_pass=m.H1_wc&&m.H2_q0&&m.H3_majority&&m.H4_weak;
  end
  res.(methods{q})=m;
 end
 out.(sp{1})=res;
end
uq_writejson(fullfile(root,'analysis','gates.json'),out);disp(jsonencode(out.val.QScanRefit));disp(jsonencode(out.ind.QScanRefit));
end
