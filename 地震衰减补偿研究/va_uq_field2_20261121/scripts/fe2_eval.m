function fe2_eval(tag)
% Evaluation with the stage-4 fld_eval.m, unchanged. fld_eval derives its input and output folders from
% the location of fld_setup.m (P.root = folder above scripts/), so both files are copied byte for byte
% (checked) into cache/evalroot/scripts/, with byte-identical copies of the stage-4 inputs it reads
% (cache/fld: meta.json, D.bin, geometry, vp.bin, inv.mat, well_F-5/F-14.bin) and this run's image file.
% fld_eval is run from that folder (the current folder precedes the path; checked with which), so
% nothing is written into uq_rtm_20261117. Outputs: analysis/fe2_eval<tag>.json (= fld_eval's
% stage4_field_fe2<tag>.json) and its two figures figures/fe2_eval_{images,metrics}<tag>.png (stage-4
% method panels).
if nargin<1,tag='';end
[C,~]=fe2_cfg();ev=fullfile(C.cache,'evalroot');fl=fullfile(ev,'cache','fld');
for d={fullfile(ev,'scripts'),fullfile(ev,'analysis'),fullfile(ev,'figures'),fl},if ~isfolder(d{1}),mkdir(d{1});end,end
for n={'fld_setup.m','fld_eval.m'}
 s=fullfile(C.s4,'scripts',n{1});t=fullfile(ev,'scripts',n{1});copyfile(s,t);
 if ~isequal(rb(s),rb(t)),error('copy of %s differs',n{1});end
end
for n={'meta.json','D.bin','shot_of.bin','rec_ix.bin','rec_iz.bin','src_ix.bin','src_iz.bin','vp.bin','inv.mat','well_F-5.bin','well_F-14.bin'}
 s=fullfile(C.s4cache,n{1});t=fullfile(fl,n{1});
 if ~isfile(t)||nb(t)~=nb(s),copyfile(s,t);end
 if ~isequal(rb(s),rb(t)),error('copy of %s differs',n{1});end
end
img=['images_fe2' tag '.mat'];copyfile(fullfile(C.cache,img),fullfile(fl,img));
old=cd(fullfile(ev,'scripts'));c=onCleanup(@() cd(old));
w1=which('fld_eval');w2=which('fld_setup');
if ~startsWith(w1,ev)||~startsWith(w2,ev),error('staged copies not first: %s | %s',w1,w2);end
fprintf('fld_eval: %s\nfld_setup: %s (byte-identical copies of uq_rtm_20261117/scripts)\n',w1,w2);
sfx=['_fe2' tag];fld_eval(img,sfx);
copyfile(fullfile(ev,'analysis',['stage4_field' sfx '.json']),fullfile(C.root,'analysis',['fe2_eval' tag '.json']));
copyfile(fullfile(ev,'figures',['field_images' sfx '.png']),fullfile(C.root,'figures',['fe2_eval_images' tag '.png']));
copyfile(fullfile(ev,'figures',['field_metrics' sfx '.png']),fullfile(C.root,'figures',['fe2_eval_metrics' tag '.png']));
delete(fullfile(fl,img));
fprintf('fe2_eval%s written %s\n',tag,datestr(now));
end

function b=rb(f),fid=fopen(f);b=fread(fid,Inf,'*uint8');fclose(fid);end

function n=nb(f),d=dir(f);n=d.bytes;end
