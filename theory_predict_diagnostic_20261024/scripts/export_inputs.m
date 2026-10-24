function export_inputs()
% Structural separation only: historical labels are known. Prediction process
% receives whitelisted observations and archived algorithm outputs, no targets.
root=fileparts(fileparts(mfilename('fullpath')));repo=fileparts(root);
src=fullfile(repo,'bae_q_20261024');addpath(fullfile(src,'core'),fullfile(src,'source'));
for model={'primary','complex','F1','F2'}
 name=model{1};calibration=[];
 if any(strcmp(name,{'primary','complex'}))
  raw=jsondecode(fileread(fullfile(repo,'soft_refit_20261023','inputs',[name '.json'])));
  ref=jsondecode(fileread(fullfile(repo,'soft_refit_20261023','results',name,'reference.json')));
  if strcmp(name,'primary'),cm=raw.config;else,cm=uq_inverse_config(raw.public,raw.estimate);end
  cm.Q=reshape(ref.reference_Q,1,[]);obs=ref.datasets;split='dev';
  if strcmp(name,'complex')
   c=jsondecode(fileread(fullfile(repo,'complex_fieldlike_20261023','core_p1','inputs','complex_input.json')));
   calibration=struct('D',c.calibrationD,'public',c.public);
  end
 else
  in=load(fullfile(src,'inputs','independent',[name '.mat']));cm=in.config;obs=in.obs;split='ind';
 end
 cm.t=cm.t(:);cm.Q=reshape(cm.Q,1,[]);cm.vp=reshape(cm.vp,1,[]);cm.h=reshape(cm.h,1,[]);
 for k=1:numel(obs)
  condition=obs(k).name;D=obs(k).D;sigma=reshape(obs(k).sigma,1,[]);
  b=load(fullfile(src,'results',split,sprintf('%s_Q+0_%s.mat',name,condition)));
  q=load(fullfile(src,'results',[split '_qscan'],sprintf('%s_Q+0_%s.mat',name,condition)));
  FullCalRefit=double(b.g.FullCalRefit);QScanRefit=double(q.g.QScanRefit);scaninfo=q.info;
  errors=[-50 -30 -10 0 10 30 50];ensemble=zeros(size(D,1),size(D,2),numel(errors));
  for e=1:numel(errors)
   z=load(fullfile(src,'results',split,sprintf('%s_Q%+d_%s.mat',name,errors(e),condition)),'g');
   ensemble(:,:,e)=double(z.g.FullCalRefit);
  end
  qspread=std(ensemble,0,3);trace_ids=unique([1:3:99,23,78]);
  save(fullfile(root,'inputs',[name '_' condition '.mat']),'name','condition','cm','D','sigma','calibration', ...
   'FullCalRefit','QScanRefit','scaninfo','qspread','trace_ids','-v7');
  fprintf('Exported solver-only %s %s\n',name,condition);
 end
end
end
