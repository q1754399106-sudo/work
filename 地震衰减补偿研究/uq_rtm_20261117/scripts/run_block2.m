function run_block2(name,B)
% Step 2 of PLAN_BLOCKFLIP_ZH.md (addendum PLAN_BLOCKFLIP_STEP2_ZH.md): imaging with the 1/snr of block length B
% (1 or B* = 300) from cache/invblk_<name>.mat. Datasets m30, m40 (data_<name>.mat) and l30
% (data_<name>_l30.mat), E1 at m = 0.5 0.7 1 1.4 2, methods of migrate2 (LP cut-off, W0, W1 follow the
% estimate). True migration velocity. Output cache/images2_<name>_bf<B>.mat, analysis/stage2_<name>_bf<B>.json.
% The E1 travel-time fields are the stage-2 ones (itag '').
P=rtm_setup();C=s2_cfg(name,1);m5=[.5 .7 1 1.4 2];
C.sets={'m30',{'E1',m5};'m40',{'E1',m5};'l30',{'E1',m5}};
C.methods={'PH','IQ30','IQ70','LP','W0','UQW0','W1','UQW1'};C.fixed_methods=true;
C.tag=sprintf('_bf%d',B);C.itag='';
C.dfiles={['data_' name '.mat'],['data_' name '_l30.mat']};
if strncmp(name,'C',1),C.dfiles=C.dfiles(1);end           % C models: l30 is in data_<name>.mat
C.invfile=sprintf('inv_%s_bf%d.mat',name,B);fi=fullfile(P.cache,C.invfile);
if ~isfile(fi)
 I=load(fullfile(P.cache,['invblk_' name '.mat']));S=struct();
 for s={'m30','m40','l30'},S.(s{1})=I.(sprintf('%s_B%d',s{1},B));end
 if B==1&&isfile(fullfile(P.cache,['inv_' name '.mat']))  % reproduction: B = 1 uses the signs of est_noise2
  O=load(fullfile(P.cache,['inv_' name '.mat']),'m30','m40');
  for s={'m30','m40'},fprintf('%s %s B1 vs inv_%s: max rel diff %.3g\n',name,s{1},name, ...
    max(abs(double(S.(s{1})(:))-double(O.(s{1})(:)))./max(double(O.(s{1})(:)),realmin)));end
 end
 save(fi,'-struct','S','-v7.3');
end
if ~isfile(fullfile(P.cache,['images2_' name C.tag '.mat'])),migrate2(C);end
score_stage2(name,1,C);fprintf('run_block2 %s B %d done\n',name,B);
end
