function run_irls(worker)
% IRLS worker: claims ../results/irls_<job>.claim, resumable. Same data as run_lasso.py (../cache).
addpath('../../qscan_orig_solver_20261106/scripts/orig_wang2021'); addpath('../../qscan_orig_solver_20261106/scripts');
L=orig_lib(); out='../results/'; ts=mean(L.windows.W3);
conds={'known','W3_0','W3_1','F0.8','F1.25'};
for s=1201:1206
  for Qt=[50 70 100 140]
    for c=1:numel(conds)
      name=sprintf('irls_Q%d_s%d__%s',Qt,s,conds{c});
      if exist([out name '.json'],'file')||exist([out name '.claim'],'file'), continue; end
      fid=fopen([out name '.claim'],'w'); fprintf(fid,'%d',worker); fclose(fid);
      m=load(sprintf('../cache/d_Q%d_s%d.mat',Qt,s));
      switch conds{c}
        case 'known', fw=L.f0;
        case 'W3_0', fw=fw_anchor(m.fp,ts,0);
        case 'W3_1', fw=fw_anchor(m.fp,ts,0.01);
        case 'F0.8', fw=0.8*L.f0;
        case 'F1.25', fw=1.25*L.f0;
      end
      t0=tic; S=scan_q2(m.dd,L,fw);
      R=struct('solver','irls','Qt',Qt,'seed',s,'cond',conds{c},'fp',m.fp,'fw',fw,'v1',S.v1,'cross',S.cross,'vd',S.vd,'v2',S.v2,'E1',S.E1,'CD',S.CD,'C2',S.C2,'vgrid',L.vgrid,'seconds',toc(t0));
      fid=fopen([out name '.json'],'w'); fprintf(fid,'%s',jsonencode(R)); fclose(fid); delete([out name '.claim']);
      fprintf('%s fw %.1f v1 %.4f (%d) vd %.4f v2 %.4f %.0fs\n',name,fw,S.v1,S.cross,S.vd,S.v2,R.seconds);
    end
  end
end
end
