% Development: dominant frequency of each wavelet window against the number of traces and the noise (seeds 901-903),
% with the value expected from the Ricker source and a(t) = t: 2/f - 2f/f0^2 - pi*t/Q = 0.
addpath('orig_wang2021'); L0=orig_lib();
[~,~,A50]=qcdn_core(zeros(L0.n1,1),L0.dt,L0.f0,1,1,50); [~,~,A100]=qcdn_core(zeros(L0.n1,1),L0.dt,L0.f0,1,1,100);
for Qt=[50 100]
  if Qt==50, A=A50; else, A=A100; end
  for ntr=[12 48 192]
    for nz=[0.01 0.05]
      L=L0; L.ntr=ntr; L.noise=nz; out=[];
      for w=fieldnames(L.windows)'
        win=L.windows.(w{1}); ts=mean(win); g=@(f) 2./f-2*f/L.f0^2-pi*ts/Qt; fe=fzero(g,[5 60]); fs=[];
        for s=901:903
          rng(s,'twister'); r=zeros(L.n1,ntr);
          for k=1:ntr, idx=find(rand(L.n1,1)<0.03); idx=idx(idx>10&idx<L.n1-10); r(idx,k)=randn(numel(idx),1); end
          x=A*r; d=x+nz*sqrt(mean(x(:).^2))*randn(size(x)); fs(end+1)=peak_freq(d,L,win); %#ok<AGROW>
        end
        out=[out sprintf(' %s exp %.1f got %s;',w{1},fe,mat2str(fs))]; %#ok<AGROW>
      end
      fprintf('Qt %d ntr %d noise %.2f:%s\n',Qt,ntr,nz,out);
    end
  end
end
