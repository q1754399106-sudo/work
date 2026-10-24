function m=score_gather(A,T,cm)
% Truth-side metrics for one gather (evaluation only).
tg=T.target;L=sqrt(movmean(tg.^2,21,1));quiet=L<.02*max(L(:));
m=struct('overall',norm(A-tg,'fro')/norm(tg,'fro'),'leak',norm(A(quiet))/norm(tg(:)), ...
 'weak',NaN,'fault',NaN,'pinch',NaN,'low_energy',NaN);
if isfield(T,'events'),m.weak=score_weak(A,T.reflectivity,T.events,T.weak,cm.fpeak,cm.dt);end
f=fieldnames(T.masks);
for i=1:numel(f),M=T.masks.(f{i});m.(f{i})=norm(A(M)-tg(M))/norm(tg(M));end
end
