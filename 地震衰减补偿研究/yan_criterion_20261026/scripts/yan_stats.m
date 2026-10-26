function yan_stats(gen,varargin)
% Per gather (PLAN_ZH.md): rebuild the noisy gather from the sweep seed (sigma must match the sweep
% bit for bit), then for each of the 13 QScan grid factors run the frozen lasso exactly as
% qscan_select.m and record the QScan objective J and the Yan 2023 statistics per trace: raw
% squared residual and sum |r|^(1/2) of the physical reflectivity r = X ./ scale. argmin J must equal
% the archived QScan choice. gen = 'G' or 'H'; further arguments are model numbers (as text).
% Writes cache/yan/<model>_<scene>_r<r>.mat. Truth is not read.
maxNumCompThreads(4);root=fileparts(fileparts(mfilename('fullpath')));in=fullfile(root,'cache','inputs');
addpath(fullfile(in,'bae_q_20261024','source'));
if gen=='G',d='qscan_offset_20261025';seed0=1300000000;else,d='fd2d_20261025';seed0=1400000000;end
addpath(fullfile(in,d,'scripts'));
sweep=fullfile(fileparts(root),d,'cache','sweep');out=fullfile(root,'cache','yan');if ~isfolder(out),mkdir(out);end
tags={'Q100_14dB','Q100_08dB','Q070_14dB'};snr=[14 8 14];qsc=[1 1 .7];lev=[1 2 1];
grid=exp(linspace(log(.6),log(2.2),13));started=tic;
for mm=varargin
 m=str2double(mm{1});name=sprintf('%s%d',gen,m);
 M=load(fullfile(in,d,'models',[name '.mat']));cm=M.config;cm.t=cm.t(:);D0=M.clean;
 for s=1:3
  for r=0:6
   f=fullfile(out,sprintf('%s_%s_r%d.mat',name,tags{s},r));if isfile(f),continue;end
   S=load(fullfile(sweep,[name '_' tags{s}],sprintf('r%03d.mat',r)),'sig','seed','info');
   if r==0,D=D0;else
    assert(S.seed==seed0+10000000*m+100000*lev(s)+r);D=D0+noise_replica(D0,snr(s),S.seed);
    assert(isequal(estimate_sigmas(D,cm.t),S.sig),'sigma mismatch %s %s r%d',name,tags{s},r);
   end
   sigma=reshape(S.sig,1,[]);[n,nx]=size(D);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;
   Q=reshape(cm.Q,1,[])*qsc(s);J=zeros(1,13);eps_=zeros(13,nx);s12=zeros(13,nx);
   for q=1:13
    c=cm;c.Q=Q*grid(q);dd=prune_dictionary(c,cm.solver_offset);rhs=dd.B'*D(:,live);
    X=refit_lasso(dd.Gram,rhs,lambda(live),dd.L);res=D(:,live)-dd.B*X;
    J(q)=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);
    eps_(q,live)=sum(res.^2,1);s12(q,live)=sum(sqrt(abs(X./dd.scale(:))),1);
   end
   [~,bestJ]=min(J);assert(bestJ==S.info.best,'QScan choice not reproduced %s %s r%d',name,tags{s},r);
   save([f '.part'],'J','eps_','s12','bestJ','live','-v7');movefile([f '.part'],f);
   fprintf('%s %s r%d QScan best %d reproduced (%.0fs)\n',name,tags{s},r,bestJ,toc(started));
  end
 end
end
fprintf('DONE %s %s\n',gen,strjoin(varargin,' '));
end
