function run_orig(worker)
% Worker of PLAN_ZH.md: claims jobs via ../results/<job>.claim, resumable; writes ../results/<job>.json.
addpath('orig_wang2021'); L=orig_lib(); out='../results/';
Qts=[50 70 100 140]; seeds=1101:1106; VA=[0 0.01]; wn={'W1','W2','W3'};
for s=seeds
  for Qt=Qts
    d=[];
    jobs={}; jobs{end+1}={'known',NaN};
    for w=1:3, for va=VA, jobs{end+1}={wn{w},va}; end, end %#ok<AGROW>
    for j=1:numel(jobs)
      cond=jobs{j}{1}; va=jobs{j}{2};
      if strcmp(cond,'known'), name=sprintf('Q%d_s%d__known',Qt,s); else, name=sprintf('Q%d_s%d__%s_va%.2f',Qt,s,cond,va); end
      if exist([out name '.json'],'file'), continue; end
      [fid,~]=fopen([out name '.claim'],'r'); if fid>0, fclose(fid); continue; end
      fid=fopen([out name '.claim'],'w'); fprintf(fid,'%d',worker); fclose(fid);
      if isempty(d), d=gen_data(L,Qt,s); end
      if strcmp(cond,'known'), fp=NaN; fw=L.f0; ts=NaN;
      else, win=L.windows.(cond); ts=mean(win); fp=peak_freq(d,L,win); fw=fw_anchor(fp,ts,va); end
      t0=tic; S=scan_q(d,L,fw);
      R=struct('Qt',Qt,'u',1/Qt,'seed',s,'cond',cond,'va',va,'ts',ts,'fp',fp,'fw',fw,'v1',S.v1,'v2',S.v2,'E1',S.E1,'C2',S.C2,'vgrid',L.vgrid,'seconds',toc(t0));
      fid=fopen([out name '.json'],'w'); fprintf(fid,'%s',jsonencode(R)); fclose(fid); delete([out name '.claim']);
      fprintf('%s fp %g fw %.1f v1 %.4f v2 %.4f %.0fs\n',name,fp,fw,S.v1,S.v2,R.seconds);
    end
  end
end
end
