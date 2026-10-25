function verify_prepared_replay(cachepath,resultpath,dest)
% Replays one observed-data trace from a binary prepared cache; no target read.
c=load(cachepath);r=jsondecode(fileread(resultpath));j=1;
out=qtrust_refit(c.D(:,j),c.d,c.sigma(j),c.initial.reflectivity(:,j),c.initial.support{j},c.dB,c.Cq,.5,true);
expected=r.Trust05(:,j);difference=norm(out.gather-expected)/max(norm(expected),realmin);
check=struct('passed',difference<1e-10,'relative_difference',difference,'trace',j,'kappa',.5,'scope','One-trace independent-process replay of selected development estimator from saved prepared cache; not independent geological validation');
uq_writejson(dest,check);assert(check.passed);
end
