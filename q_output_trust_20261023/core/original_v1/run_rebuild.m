function run_rebuild()
% Entry point: run_rebuild. Consult VALIDATION.md for execution status.
c=config_model(); if ~exist(c.output,'dir'), mkdir(c.output); end
rng(c.seed); n=c.nt; nx=numel(c.offset); I=eye(n);
[a0,~,~,~]=ray_absorption(c.t,0,c); G0=attenuation_matrix(a0,c);
B0=ridge_operator(G0,c.lambda_fixed);
names={'Truth','Observed','Step1','FixedSplit','AdaptiveUnweighted','CovarianceSplit','Direct'};
D=zeros(n,nx,numel(names)); rows=[]; diagrows=[]; qrows=[];
for ix=1:nx
    x=c.offset(ix); fprintf('Offset %g m (%d/%d)\n',x,ix,nx);
    [~,ax,as]=ray_absorption(c.t,x,c);
    Gx=attenuation_matrix(ax,c); Gt=attenuation_matrix(a0+ax,c);
    s=make_reference(x,c); y0=G0*s; clean=Gt*s;
    sigma=norm(clean)/sqrt(n)/10^(c.snr_db/20);
    noise=sigma*randn(n,1); d=clean+noise;
    Bx=ridge_operator(Gx,c.lambda_x); y=Bx*d;
    fixed=B0*y;
    % Covariance weighting uses synthetic known sigma: an oracle-noise benchmark.
    Cunit=Bx*Bx'; L=chol(Cunit+c.cov_floor*I,'lower');
    target=sigma*sqrt(n); % near-exact with a small covariance floor
    [weighted,lw]=discrepancy_solve(L\G0,L\y,c.lambda_grid,target);
    [adaptive,la]=discrepancy_solve(G0,y,c.lambda_grid,sigma*norm(Bx,'fro'));
    [direct,ld]=discrepancy_solve(Gt,d,c.lambda_grid,target);
    out=[s,d,y,fixed,adaptive,weighted,direct]; D(:,ix,:)=reshape(out,n,1,[]);
    for k=2:7
        ref=s; if k==3, ref=y0; end
        err=norm(out(:,k)-ref)/norm(ref);
        rows=[rows; x,k,err,20*log10(norm(ref)/max(norm(out(:,k)-ref),eps)),lw,la,ld]; %#ok<AGROW>
    end
    e1=y-y0;
    identity=norm((fixed-s)-((B0*G0-I)*s+B0*e1))/max(norm(s),eps);
    idx=c.t>=0.25 & c.t<=1.05;
    smallerr=norm(as(idx)-ax(idx))/max(norm(ax(idx)),eps);
    mismatch=norm(Gt-Gx*G0,'fro')/norm(Gt,'fro');
    noisegain=norm(B0*Bx*noise)/max(norm(noise),eps);
    diagrows=[diagrows;x,smallerr,mismatch,norm(e1)/norm(y0),noisegain,identity]; %#ok<AGROW>
    % Q-mismatch subset: same observed trace, same physical ray geometry.
    if ix==1 || ix==nx
        for qs=c.q_scales
            cm=c; cm.Q=c.Q*qs;
            [aa,xx]=ray_absorption(c.t,x,cm);
            A=attenuation_matrix(aa,cm); X=attenuation_matrix(xx,cm);
            BB=ridge_operator(X,c.lambda_x); yy=BB*d;
            LL=chol(BB*BB'+c.cov_floor*I,'lower');
            [ss,ll]=discrepancy_solve(LL\A,LL\yy,c.lambda_grid,target);
            qrows=[qrows;x,qs,norm(ss-s)/norm(s),ll]; %#ok<AGROW>
        end
    end
end
T=array2table(rows,'VariableNames',{'Offset_m','MethodIndex','RelativeError','OutputSNR_dB','LambdaCov','LambdaAdaptive','LambdaDirect'});
T.Method=string(names(T.MethodIndex))'; writetable(T,fullfile(c.output,'metrics.csv'));
J=array2table(diagrows,'VariableNames',{'Offset_m','SmallAngleRelativeError','MatrixCompositionError','Step1RelativeError','RealizedNoiseGain','ErrorIdentityResidual'});
writetable(J,fullfile(c.output,'diagnostics.csv'));
E=event_metrics(D,names,c); writetable(E,fullfile(c.output,'event_metrics.csv'));
writetable(array2table(qrows,'VariableNames',{'Offset_m','QScale','RelativeError','Lambda'}),fullfile(c.output,'q_sensitivity.csv'));
save(fullfile(c.output,'experiment.mat'),'c','D','names','T','J','qrows');
fig=figure('Visible','off','Position',[100 100 1400 750]); tiledlayout(2,3);
sel=[1 2 3 4 6 7]; lim=max(abs(D(:,:,1)),[],'all');
for k=sel
 nexttile; imagesc(c.offset,c.t,D(:,:,k),[-lim lim]); title(names{k}); xlabel('Offset (m)'); ylabel('Time (s)'); ylim([0.2 1.1]); colorbar;
end
colormap(gray); exportgraphics(fig,fullfile(c.output,'gathers.png'),'Resolution',150); close(fig);
fig=figure('Visible','off','Position',[100 100 1000 450]); tiledlayout(1,2);
nexttile; hold on;
for k=[2 4 5 6 7], z=T.MethodIndex==k; plot(T.Offset_m(z),T.RelativeError(z),'-o','DisplayName',names{k}); end
legend('Location','best'); xlabel('Offset (m)'); ylabel('Relative waveform error'); grid on;
nexttile; plot(J.Offset_m,J.MatrixCompositionError,'-o'); xlabel('Offset (m)'); ylabel('||Gt-GxG0||F / ||Gt||F'); grid on;
exportgraphics(fig,fullfile(c.output,'errors.png'),'Resolution',150); close(fig);
run_checks(c,G0,J);
fprintf('Completed. Outputs: %s\n',c.output);
disp(groupsummary(T,'Method','mean','RelativeError'));
end

function run_checks(c,G0,J)
n=c.nt; I=eye(n); Z=attenuation_matrix(zeros(n,1),c);
assert(norm(Z-I,'fro')/sqrt(n)<1e-11,'Zero absorption must be identity');
[~,ax]=ray_absorption(c.t,0,c); assert(max(abs(ax))==0);
[aa,xx]=ray_absorption(c.event_t,1000,c);
assert(all(aa>0) && all(xx>=0),'Ray absorption must be nonnegative');
assert(max(J.ErrorIdentityResidual)<1e-10,'Error decomposition identity failed');
u=sin((1:n)'); v=cos((1:n)');
assert(abs((G0*u)'*v-u'*(G0'*v))<1e-10,'Adjoint check failed');
assert(all(isfinite(J.MatrixCompositionError)));
fid=fopen(fullfile(c.output,'checks.txt'),'w');
fprintf(fid,'PASS: zero-Q-loss identity, zero-offset relative loss, nonnegative path loss, adjoint, step-error identity.\n'); fclose(fid);
end
