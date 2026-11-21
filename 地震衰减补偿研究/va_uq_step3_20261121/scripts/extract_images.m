function extract_images(name)
% Light extraction (no migration, no LU) for va_uq_step3_20261121: copies the interior images of one
% development model from the existing caches into a scipy-readable v7 file cache/img_<name>.mat.
%   D2, D4: va_uq_step1_20261121/cache/images_va_<name>.mat (v_band migration; R0 and the step-1
%           methods on m = 0.5 0.7 1 1.4 2), reference = elastic data migrated with v100
%           (uq_rtm_20261117/cache/images2_<name>.mat, ref.R0), model v, Q, v_band, anchor z_a(x).
%   D1, D3: uq_rtm_20261117/cache/images2_<name>.mat (stage 2, v100 migration; m30 R0, PH, UQW1 on the
%           E1 grid 0.5 ... 2), same reference; model v, Q, vmig; anchor = first macro interface.
% Only the files are read; nothing in other directories is written. PLAN_ZH.md section 5.
here=fileparts(mfilename('fullpath'));root=fileparts(here);base=fileparts(root);
s1=fullfile(base,'va_uq_step1_20261121','scripts');src=fullfile(base,'uq_rtm_20261117','scripts');
addpath(src);addpath(s1);
P=rtm_setup();M=make_model(name,P);iz=P.in{1};ix=P.in{2};
O=load(fullfile(P.cache,['images2_' name '.mat']),'ref');
S=struct('name',name,'dx',P.dx,'z',P.z(iz),'x',P.x(ix),'v',M.v(iz,ix),'Q',M.Q(iz,ix),'vmig',M.vmig(iz,ix), ...
 'ref',reshape(double(O.ref.R0),P.nzi,P.nxi));clear O
switch name
 case {'D2','D4'}
  [C,~]=va_cfg(name);vb=va_vband(M,P,C);S.vb=vb(iz,ix);S.za=M.spec.iface{C.anchor_iface}(P.x(ix));
  I=load(fullfile(base,'va_uq_step1_20261121','cache',['images_va_' name '.mat']));
  S.m=I.meta.m;S.methods=I.meta.methods;S.R0=reshape(double(I.R0),P.nzi,P.nxi);
  for q=I.meta.methods,S.(q{1})=reshape(double(I.(q{1})),P.nzi,P.nxi,[]);end
 case {'D1','D3'}
  S.vb=S.vmig;S.za=M.spec.iface{1}(P.x(ix));
  I=load(fullfile(P.cache,['images2_' name '.mat']),'m30','meta');rows=I.meta.combo.m30;I=I.m30;
  e1=find(strcmp(rows(:,1),'E1'))';S.m=cell2mat(rows(e1,2))';S.methods={'PH','UQW1'};
  S.R0=reshape(double(I.R0),P.nzi,P.nxi);
  for q=S.methods,S.(q{1})=reshape(double(I.(q{1})(:,e1)),P.nzi,P.nxi,[]);end
 otherwise
  error('development models only (D1-D4)');
end
save(fullfile(root,'cache',['img_' name '.mat']),'-struct','S','-v7');
fprintf('extract %s: %s\n',name,strjoin(fieldnames(S)',' '));
end
