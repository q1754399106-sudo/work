function run_baselines(split)
% Experiment 2 literature baselines on the same scenes, sigma and applied Q.
% split='dev' uses soft_refit development scenes; 'val' the frozen V1-V3 inputs.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
opts=jsondecode(fileread(fullfile(root,'SPEC.json')));gains=[10 20 30];started=tic;
outdir=fullfile(root,'results',[split '_baselines']);if ~isfolder(outdir),mkdir(outdir);end
scenes=list_scenes(root,split);
for k=1:numel(scenes)
 sc=scenes(k);
 for errorpct=reshape(opts.errors_percent,1,[])
  file=fullfile(outdir,sprintf('%s_Q%+d_%s.mat',sc.model,errorpct,sc.name));if isfile(file),continue;end
  offset=450;if isfield(sc.cm,'solver_offset'),offset=sc.cm.solver_offset;end
  used=sc.cm;used.Q=reshape(sc.referenceQ,1,[])*(1+errorpct/100);d=prune_dictionary(used,offset);
  g=struct();sbl=zeros(size(sc.D));conv=true;
  for first=1:11:size(sc.D,2)
   ids=first:min(first+10,size(sc.D,2));o=sbl_tvd(sc.D(:,ids),d,sc.sigma(ids));sbl(:,ids)=o.gather;conv=conv&&o.converged;
  end
  g.SBL_TVD=single(sbl);
  for G=gains,g.(sprintf('WangIQF_G%02d',G))=single(wang_inverse_q(sc.D,used,offset,G));end
  info=struct('used_Q',used.Q,'sbl_converged',conv); %#ok<NASGU>
  save(file,'g','info','-v7');
  fprintf('%s Q%+d %s sbl_converged=%d elapsed=%.1fs\n',sc.model,errorpct,sc.name,conv,toc(started));
 end
end
fprintf('DONE %s %.1fs\n',split,toc(started));
end
