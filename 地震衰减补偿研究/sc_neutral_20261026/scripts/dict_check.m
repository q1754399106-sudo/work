function dict_check()
% Pre-run input check: exports the frozen MATLAB dictionary (prune_dictionary, as used by the
% Standard pipeline of fd2d) for H1-H6 at Q0 and Q0x0.7, so that scripts/check_dictionary.py can
% compare it with the Python dictionary used by SC. No method output is computed.
root=fileparts(fileparts(mfilename('fullpath')));in=fullfile(root,'cache','inputs');
addpath(fullfile(in,'bae_q_20261024','source'));out=fullfile(root,'cache','dictcheck');if ~isfolder(out),mkdir(out);end
for m=1:6
 M=load(fullfile(in,'fd2d_20261025','models',sprintf('H%d.mat',m)));cm=M.config;cm.t=cm.t(:);
 for qs=[1 .7]
  c=cm;c.Q=reshape(cm.Q,1,[])*qs;d=prune_dictionary(c,cm.solver_offset);B=d.B;W=d.W;scale=d.scale;
  save(fullfile(out,sprintf('H%d_q%03d.mat',m,round(100*qs))),'B','W','scale','-v7');
 end
end
fprintf('DONE dict_check\n');
end
