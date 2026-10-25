function run_crit(first,last)
% Worker (PLAN_ZH.md): scenes first..last, replicates 1-6. Regenerates each observation from the
% seed stored in the sweep cache, then evaluates A0-A3 (qscan_criterion_20261024) and C1-C4
% (qcriteria_20261024) on the 13 QScan grid factors, all traces. A0 must reproduce the QScan
% choice stored in the cache. Truth is not read.
P=cc_setup();maxNumCompThreads(P.threads);out=fullfile(P.root,'cache','crit');if ~isfolder(out),mkdir(out);end
ng=numel(P.grid);started=tic;
for s=first:last
 sc=P.scenes(s);[cm,D0]=load_model(P,sc);used=cm;used.Q=reshape(cm.Q,1,[])*sc.q_scale;offset=cm.solver_offset;
 for r=1:P.R
  file=fullfile(out,sprintf('%s_r%03d.mat',sc.name,r));if isfile(file),continue;end
  t0=tic;L=load(fullfile(sc.sweep,sprintf('r%03d.mat',r)),'info','seed','sig');
  D=D0+noise_replica(D0,sc.snr_db,L.seed);sig=reshape(L.sig,1,[]);
  st=qscan_stats(D,sig,used,offset,P.grid);st=rmfield(st,'gathers');
  cv=zeros(4,ng);
  for q=1:ng
   c=used;c.Q=used.Q*P.grid(q);d=prune_dictionary(c,offset);cv(:,q)=criteria(D,d.G*d.W,d,sig);
  end
  sel=zeros(1,8);score=nan(8,ng);
  for k=1:4,[sel(k),score(k,:)]=select_factor(st,P.criteria{k});end
  for k=2:4,score(3+k,:)=cv(k,:);[~,sel(3+k)]=min(cv(k,:));end
  c4=spectral_factor(D,used,offset,P.grid);sel(8)=find(P.grid==c4,1);
  assert(sel(1)==L.info.best,'A0 does not reproduce the cached QScan choice (%s r=%d).',sc.name,r);
  assert(all(abs(cv(1,:)'-score(1,:)')<=1e-9*abs(score(1,:)')),'C0 and A0 objectives differ (%s r=%d).',sc.name,r);
  info=struct('criteria',{P.criteria},'grid',P.grid,'optimization_flags',st.optimization_flags,'seconds',toc(t0));
  save([file '.part'],'sel','score','c4','info','-v7');movefile([file '.part'],file);
  fprintf('%s r=%d sel %s flags %d (%.0fs, total %.0fs)\n',sc.name,r,mat2str(sel),st.optimization_flags,info.seconds,toc(started));
 end
end
fprintf('DONE scenes %d-%d %.0fs\n',first,last,toc(started));
end
