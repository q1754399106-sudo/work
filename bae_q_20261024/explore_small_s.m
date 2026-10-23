function explore_small_s(which,s_values)
% POST-HOC exploration (decided after seeing development results; see RESULTS):
% BAE with smaller prior widths on development scenes only.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
opts=jsondecode(fileread(fullfile(root,'SPEC.json')));soft=fullfile(fileparts(root),'soft_refit_20261023');[cm,ref]=load_dev(soft,which);
outdir=fullfile(root,'results','dev_explore');if ~isfolder(outdir),mkdir(outdir);end;started=tic;
for errorpct=reshape(opts.errors_percent,1,[])
 used=cm;used.Q=reshape(ref.reference_Q,1,[])*(1+errorpct/100);d=prune_dictionary(used,450);
 ops=cell(size(s_values));for q=1:numel(s_values),ops{q}=q_sample_operators(used,450,s_values(q),opts.M,opts.rho,opts.seed);end
 for k=1:numel(ref.datasets)
  data=ref.datasets(k);file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',which,errorpct,data.name));if isfile(file),continue;end
  D=data.D;sig=reshape(data.sigma,1,[]);g=struct();for q=1:numel(s_values),g.(sprintf('BAE_s%03d',round(100*s_values(q))))=zeros(size(D),'single');end
  for first=1:11:size(D,2)
   ids=first:min(first+10,size(D,2));initial=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sig(ids));
   for q=1:numel(s_values),b=bae_q_refit(D(:,ids),d,ops{q},sig(ids),initial,true);g.(sprintf('BAE_s%03d',round(100*s_values(q))))(:,ids)=single(b.gather);end
  end
  save(file,'g','-v7');fprintf('%s Q%+d %s %.1fs\n',which,errorpct,data.name,toc(started));
 end
end
fprintf('DONE explore %s\n',which);
end
