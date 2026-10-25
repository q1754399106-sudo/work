function prep_noise()
% Rebuild the noisy H1-H6 gathers of fd2d_20261025 (PLAN_ZH.md, 输入). Seeds, SNR and noise process
% are those of fd2d_20261025/scripts/fd_setup.m and run_fd.m: level 1 = 14 dB, level 2 = 8 dB,
% seed 1.4e9 + 1e7*m + 1e5*level + r, r = 1..6; r = 0 is the noise-free gather. The rebuilt gather is
% accepted only if estimate_sigmas reproduces the sigma saved in the fd2d sweep file exactly (r >= 1)
% and the saved seed matches. Writes cache/noisy/H<m>_L<level>_r<r>.mat (D, sig) and
% logs/prep_noise.json. Truth is not read.
% Four threads as in run_fd.m: with one thread the sums differ in the last bit (~1e-15 relative).
maxNumCompThreads(4);root=fileparts(fileparts(mfilename('fullpath')));in=fullfile(root,'cache','inputs');
addpath(fullfile(in,'fd2d_20261025','scripts'),fullfile(in,'bae_q_20261024','source'));
sweep=fullfile(fileparts(root),'fd2d_20261025','cache','sweep');out=fullfile(root,'cache','noisy');
if ~isfolder(out),mkdir(out);end
tags={'Q100_14dB','Q100_08dB'};snr=[14 8];checks={};
for m=1:6
 M=load(fullfile(in,'fd2d_20261025','models',sprintf('H%d.mat',m)));D0=M.clean;t=M.config.t(:);
 for level=1:2
  for r=0:6
   S=load(fullfile(sweep,sprintf('H%d_%s',m,tags{level}),sprintf('r%03d.mat',r)),'sig','seed');
   if r==0
    D=D0;seed=NaN;dsig=NaN;
   else
    seed=1400000000+10000000*m+100000*level+r;assert(S.seed==seed);
    D=D0+noise_replica(D0,snr(level),seed);dsig=max(abs(estimate_sigmas(D,t)-S.sig));assert(dsig==0);
   end
   sig=S.sig;save(fullfile(out,sprintf('H%d_L%d_r%d.mat',m,level,r)),'D','sig','-v7');
   checks(end+1,:)={sprintf('H%d',m),tags{level},r,seed,dsig}; %#ok<AGROW>
  end
 end
 % Q0x0.7 / 14 dB reuses the level-1 noise: its saved sigma must equal the Q0 / 14 dB one.
 for r=0:6
  A=load(fullfile(sweep,sprintf('H%d_Q100_14dB',m),sprintf('r%03d.mat',r)),'sig');
  B=load(fullfile(sweep,sprintf('H%d_Q070_14dB',m),sprintf('r%03d.mat',r)),'sig');assert(isequal(A.sig,B.sig));
 end
end
T=cell2table(checks,'VariableNames',{'model','scene','r','seed','max_sigma_difference'});disp(T);
fid=fopen(fullfile(root,'logs','prep_noise.json'),'w');fprintf(fid,'%s',jsonencode(table2struct(T),'PrettyPrint',true));fclose(fid);
end
