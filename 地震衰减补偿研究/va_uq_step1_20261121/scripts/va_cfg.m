function [C,P]=va_cfg(name)
% Configuration of va_uq_step1_20261121 (PLAN_ZH.md). Reuses the engine of uq_rtm_20261117 (addpath,
% unchanged); inputs (data, noise estimate, reference image) are read from its cache, everything
% computed here goes to this experiment's own cache/.
here=fileparts(mfilename('fullpath'));C.root=fileparts(here);
src=fullfile(fileparts(C.root),'uq_rtm_20261117','scripts');if ~contains(path,src),addpath(src);end
P=rtm_setup();C.src_cache=P.cache;P.cache=fullfile(C.root,'cache');   % own cache for t* fields, state, images
if ~isfolder(P.cache),mkdir(P.cache);end
C.name=name;C.tag='';
C.dataset='m30';C.m=[.5 .7 1 1.4 2];            % E1: Q_c = m Q_m
C.f_band=25;                                    % migration velocity = phase velocity at 25 Hz
C.fv=[15 35];C.fv_bar=sqrt(prod(C.fv));C.Lam=log(C.fv(2)/C.fv(1));   % VA-UQ: ln f_v uniform on [15,35] Hz
C.sig=P.sig;                                    % 0.447, on Delta t* only (VA) / on t* (UQW1)
switch name                                     % anchor: macro interface (index into M.spec.iface)
 case 'D2',C.anchor_iface=1;                    % 450 - 120 exp(-((x-1200)/500)^2): 330-450 m
 case 'D4',C.anchor_iface=2;                    % flat 330 m
 otherwise,error('only the development models D2, D4 are allowed in this experiment');
end
% methods: label -> {table, t* key}; PH25, VAfix are descriptive only (not in P1-P3)
C.methods={'PH100','UQW1_100','UQW1_25','VA','PH25','VAfix'};
C.win=struct('deep',{{[300 2100],[820 1550]}});
switch name
 case 'D2',C.win.body={[1050 1350],[850 1300]};
 case 'D4',C.win.body={[300 2100],[400 800]};
end
P.threads=12;
end
