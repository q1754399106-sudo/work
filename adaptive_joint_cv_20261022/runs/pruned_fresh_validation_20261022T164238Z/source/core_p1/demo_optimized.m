function result=demo_optimized()
root=fileparts(mfilename('fullpath'));a=load(fullfile(root,'examples','example_inputs.mat'));
result=compensate_seismic(a.D,a.offsets,a.public,a.calibrationD,a.calibrationOffsets,a.domain);
save(fullfile(root,'examples','your_demo_output.mat'),'result','-v7');
fprintf('FullCalRefit completed. Gather: %d x %d; estimated Q: %s\n',size(result.gather,1),size(result.gather,2),mat2str(result.estimated_Q.layer_Q,5));
if nargout==0
 figure('Name','Unified seismic compensation');tiledlayout(1,2);limit=max(abs([a.D(:);result.gather(:)]));
 nexttile;imagesc(1:size(a.D,2),a.public.t,a.D);clim([-limit,limit]);title('Observed');xlabel('Trace');ylabel('Time (s)');colorbar;
 nexttile;imagesc(1:size(a.D,2),a.public.t,result.gather);clim([-limit,limit]);title('Optimized');xlabel('Trace');ylabel('Time (s)');colorbar;
end
end
