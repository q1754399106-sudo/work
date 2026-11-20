function P=ou_setup()
% Constants of ou_uqrtm_20261120 (PLAN_ZH.md). Grid, acquisition, frequencies and UQ constants are those
% of uq_rtm_20261117 (rtm_setup, read through addpath). Inputs of that experiment (data, t* fields,
% noise estimates, images) are read from its cache in the main checkout and never written.
here=fileparts(mfilename('fullpath'));root=fileparts(here);
src=fullfile(fileparts(root),'uq_rtm_20261117','scripts');addpath(src);
P=rtm_setup();
P.root=root;
P.src_cache=getenv('OU_SRC_CACHE');
if isempty(P.src_cache),P.src_cache='D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache';end
P.cache=fullfile(root,'cache');                               % own outputs (not uploaded)
P.threads=3;                                                  % 5 sessions share 16 cores today
% random-field prior of the compensation error du = Q_true^-1/Q_model^-1 - 1 (PLAN_ZH.md section 3)
P.ou=struct('Lz',50,'Lx',250,'s_bg',.15,'s_flag',.6,'nseg',48,'dec',4);
end
