function fig_gathers()
% Descriptive: one shot gather of P1 (shot 26, x = 1200 m) before migration: elastic reference,
% attenuated, attenuated + noise (-30 dB), and route-C compensated (per-trace UQ-W, m = 1 and 0.7).
% Band 3-55 Hz from the stored frequency samples (inverse DFT, display only). figures/gathers_P1.png
P=rtm_setup();k=26;D=load(fullfile(P.cache,'data_P1.mat'),'ref','att','m30','C_m100','C_m070');
t=P.t;E=exp(1i*2*pi*t*P.f);off=(P.ixr-P.ixs(k))*P.dx;
pan={'ref','参照（无衰减、无噪声）';'att','有衰减（Q 模型）';'m30','有衰减 + 道上 −30 dB 噪声';
     'C_m100','逐道 UQ-W 补偿，m = 1';'C_m070','逐道 UQ-W 补偿，m = 0.7'};
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
f=figure('Visible','off','Position',[0 0 1500 560]);tl=tiledlayout(1,5,'TileSpacing','compact','Padding','compact');
G=cell(1,5);
for q=1:5,G{q}=real(E*double(D.(pan{q,1})(:,:,k)).');end
tg=t>=.2;
for q=1:5
 g=G{q};g=g./max(std(g(t>=.9&t<=1.5,:),0,1),eps);          % per-trace balance on 0.9-1.5 s, display only
 nexttile;imagesc(off,t,g,[-3 3]);colormap(gray);title(pan{q,2});xlabel('偏移距（m）');
 if q==1,ylabel('时间（s）');end
 ylim([.2 1.8]);
end
title(tl,'P1 第 26 炮（x = 1200 m）偏移前的共炮点道集，3–55 Hz；每道按 0.9–1.5 s 的能量归一，只为显示');
exportgraphics(f,fullfile(P.root,'figures','gathers_P1.png'),'Resolution',160);close(f);
% amplitude spectra of near-offset traces, 0.9-1.5 s
f=figure('Visible','off','Position',[0 0 700 380]);hold on;near=abs(off)<=300;w=t>=.9&t<=1.5;fq=(0:255)/(512*P.dt);
cl=[0 0 0;.85 .55 .2;.6 .6 .6;.8 .15 .15;.3 .55 .8];
for q=1:5
 X=abs(fft(G{q}(w,near).*hann(nnz(w)),512));s=mean(X(1:256,:),2);
 plot(fq,20*log10(s/max(s)),'Color',cl(q,:),'LineWidth',1.6,'DisplayName',pan{q,2});
end
xlim([0 70]);ylim([-50 3]);grid on;box on;xlabel('频率（Hz）');ylabel('振幅谱（dB，各自归一）');legend('Location','southwest');
title('近偏移距（|h| ≤ 300 m）0.9–1.5 s 的平均振幅谱');
exportgraphics(f,fullfile(P.root,'figures','gathers_P1_spectra.png'),'Resolution',160);close(f);
end
