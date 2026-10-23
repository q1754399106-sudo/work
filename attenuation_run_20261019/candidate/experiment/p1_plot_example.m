function p1_plot_example(example,outdir,stem)
% Only invisible, owned figures are created; existing figures are untouched.
fields=fieldnames(example.methods); include={'Frozen_Split','CombinedUncertainty_Direct','CombinedUncertainty_Split'};
fields=include(ismember(include,fields));
f=figure('Visible','off','Color','w','Position',[80 80 1450 760]);
cleanup=onCleanup(@()delete(f)); %#ok<NASGU>
tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
lim=max(abs(example.truth(:))); limit=[-lim lim];
nexttile; imagesc(1:size(example.truth,2),example.time,example.truth); caxis(limit); colorbar;
title('真值 / Truth'); xlabel('道 / Trace'); ylabel('时间 / Time (s)');
nexttile; imagesc(1:size(example.observation,2),example.time,example.observation); caxis(limit); colorbar;
title('观测 / Observation'); xlabel('道 / Trace'); ylabel('时间 / Time (s)');
for k=1:min(numel(fields),3)
 nexttile; imagesc(1:size(example.truth,2),example.time,example.methods.(fields{k})); caxis(limit); colorbar;
 title(strrep(fields{k},'_',' / '),'Interpreter','none'); xlabel('道 / Trace'); ylabel('时间 / Time (s)');
end
nexttile; trace=example.fixed_trace_indices(2); plot(example.time,example.truth(:,trace),'k','LineWidth',1.4); hold on;
labels={'真值 / Truth'};
for k=1:numel(fields), plot(example.time,example.methods.(fields{k})(:,trace)); labels{end+1}=strrep(fields{k},'_',' / '); end %#ok<AGROW>
xlabel('时间 / Time (s)'); ylabel('振幅 / Amplitude'); title(sprintf('固定第%d道 / Fixed trace %d',trace,trace));
legend(labels,'Interpreter','none','Location','best'); grid on;
sgtitle(sprintf('%s | %s | %s | cal %d | seed %d',example.Model,example.Suite,example.Scenario,example.CalibrationSeed,example.Seed),'Interpreter','none');
exportgraphics(f,fullfile(outdir,[stem '_gathers.png']),'Resolution',160);
f2=figure('Visible','off','Color','w','Position',[80 80 1200 620]);
cleanup2=onCleanup(@()delete(f2)); %#ok<NASGU>
tiledlayout(f2,1,2,'Padding','compact');
nexttile; hold on; labels={}; trace=example.fixed_trace_indices(2); dt=mean(diff(example.time));
targets=[{'truth'},fields];
for k=1:numel(targets)
 name=targets{k}; S=example.truth; if k>1, S=example.methods.(name); end
 spectrum=abs(fft(S(:,trace))); hz=(0:size(S,1)-1)/(size(S,1)*dt);
 plot(hz,spectrum); labels{end+1}=strrep(name,'_',' / '); %#ok<AGROW>
end
xlim([0 200]); xlabel('频率 / Frequency (Hz)'); ylabel('幅度谱 / Magnitude spectrum');
title('固定道频谱 / Fixed-trace spectra'); legend(labels,'Interpreter','none'); grid on;
nexttile; hold on; plot(example.time,example.truth(:,trace),'k','LineWidth',1.4);
for k=1:numel(fields), plot(example.time,example.methods.(fields{k})(:,trace)-example.truth(:,trace)); end
xlabel('时间 / Time (s)'); ylabel('误差 / Error'); title('固定道误差 / Fixed-trace error');
legend([{'真值 / Truth'},strcat(fields,' error')],'Interpreter','none'); grid on;
exportgraphics(f2,fullfile(outdir,[stem '_spectra_error.png']),'Resolution',160);
end
