function make_validation_data()
% Builds the three registered validation sets. Solver inputs and truth are
% written to separate files; the solver side never loads truth/.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
assert(isfile(fullfile(root,'analysis','SELECTED.json')),'Freeze the development selection first.');
out=fullfile(root,'inputs','validation');truth=fullfile(root,'truth');
assert(~isfolder(out)&&~isfolder(truth),'Validation data already exist; refusing to regenerate.');mkdir(out);mkdir(truth);
softroot=fullfile(fileparts(root),'soft_refit_20261023');primary=jsondecode(fileread(fullfile(softroot,'inputs','primary.json')));
cm=primary.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=cm.Q(:)';dt=cm.dt;nx=99;nt=cm.nt;
% V1: original four-layer reflectivity, new noise seeds.
cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{2026102401,2026102402});
ds=gen_layered_scene(cm,primary.reflectivity,cond);
save_layered(out,truth,'V1',cm,ds,primary.reflectivity,primary.event_indices_zero_based+1,[2 4 6 8],false(nt,nx),false(nt,nx),'original four-layer geometry, new noise');
% V2: FaultPinch four-layer variant.
cm2=cm;cm2.Q=[55 45 80 120];u=linspace(0,10,nx);
times0=[.22 .31 .45 .58 .74 .90 1.06 1.30 1.52];amps=[.09 .011 -.05 -.008 .075 .010 -.04 -.007 .03];
R=zeros(nt,nx);ev=nan(nx,9);fault=false(nt,nx);pinch=false(nt,nx);t=cm.t;
for j=1:nx
 ti=times0+.030*(j-1)/98+.006*sin(u(j)*.5+(0:8)*.8);
 if j>=50,ti(times0>=.55)=ti(times0>=.55)+.028;end
 sep=interp1([1 45 80 99],[.16 .16 .035 .035],j);ti(6)=ti(5)+sep;
 a=amps.*(1+.10*sin(u(j)*.6+(0:8)));a(6)=a(6)*min(max((80-j)/20,0),1);
 ii=round(ti/dt)+1;ev(j,:)=ii;ev(j,a==0)=NaN;
 for e=1:9,if a(e)~=0,R(ii(e),j)=R(ii(e),j)+a(e);end,end
 if j>=46&&j<=54,fault(t>=.50,j)=true;end
 if j>=60,pinch(abs(t-ti(6))<=.04,j)=true;end
end
cond2=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{2026102411,2026102412});
ds=gen_layered_scene(cm2,R,cond2);
save_layered(out,truth,'V2',cm2,ds,R,ev,[2 4 6 8],fault,pinch,'FaultPinch: dip, 28 ms fault throw at trace 50, weak reflector 6 thins and vanishes by trace 80, Q=[55 45 80 120]');
% V3: new complex geology, independent calibration with frozen layered code.
g=struct('seed',2026102421,'calibration_seed',2026102422,'vp0',1880,'q0',33,'q_gradient',.038, ...
 'fold_amp',[38 22],'fold_phase',1.1,'faults',[900 2010],'throws',[-130 120],'steps',[620 110;1000 120;1500 140], ...
 'channel',[2250 300 760],'wedge_start',300,'wedge_tip',1650,'wedge_depth',1120,'wedge_thickness',48,'gas',[600 1260], ...
 'diffractors',[900 .74 .016;2010 1.00 .013;2250 .70 .010],'dead',[31 66],'burst_traces',[12 47 88], ...
 'groundroll_x',2400,'linear_t0',.55,'linear_slope',-.00008, ...
 'observations',struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{2026102431,2026102432}, ...
  'coherent_fraction',{.18,.30},'multiple_scale',{1,1.5},'static_std_ms',{1.5,3}));
started=tic;model=gen_complex_model(g);fprintf('V3 forward model ready %.1fs\n',toc(started));
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));p=model.public;
[est,corrected,cal]=p1_prepare_candidate(model.calibrationD,p.offset,p,policy);cmfull=uq_inverse_config(corrected,est);
fprintf('V3 calibration done %.1fs: layer Q %s, f %.2f Hz\n',toc(started),mat2str(cmfull.Q,4),cmfull.fpeak);
[fault3,pinch3,lowE,deep]=complex_regions(model,g);
obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
for k=1:numel(model.datasets)
 D=model.datasets(k).D;sigma=zeros(1,nx);
 for j=1:nx,if all(D(:,j)==0),sigma(j)=realmin;else,ni=estimate_trace_snr(D(:,j),model.t);sigma(j)=ni.sigma;end,end
 obs(k)=struct('name',model.datasets(k).name,'D',D,'sigma',sigma);
 tr(k)=struct('name',model.datasets(k).name,'target',model.datasets(k).target);
end
config=cmfull;calibration=struct('estimate',est,'diagnostics',cal);
save(fullfile(out,'V3.mat'),'config','obs','calibration','-v7');
mt=model.truth;truthQ=mt.calibration_Q;regions=struct('fault',fault3,'pinch',pinch3,'low_energy',lowE,'deep',deep);
save(fullfile(truth,'V3_truth.mat'),'tr','regions','mt','g','truthQ','-v7');
uq_writejson(fullfile(out,'V3_summary.json'),struct('calibrated_Q',cmfull.Q,'calibrated_frequency',cmfull.fpeak, ...
 'coarse_true_harmonic_Q',truthQ,'q_boundary',isfield(est,'q_boundary')&&any(est.q_boundary),'description','new complex geology (MATLAB port), seeds 2026102421/22/31/32'));
fprintf('Validation data written %.1fs\n',toc(started));
end

function save_layered(out,truth,name,cm,ds,R,events,weak,fault,pinch,description)
config=cm;obs=rmfield(ds,'target');tr=struct('name',{ds.name},'target',{ds.target});
save(fullfile(out,[name '.mat']),'config','obs','-v7');
regions=struct('fault',fault,'pinch',pinch);reflectivity=R; %#ok<NASGU>
save(fullfile(truth,[name '_truth.mat']),'tr','reflectivity','events','weak','regions','description','-v7');
end

function [fault,pinch,lowE,deep]=complex_regions(model,g)
% Truth-derived scoring masks (evaluation only).
t=model.t;x=model.x;[nt,nx]=size(model.datasets(1).target);mt=model.truth;
fault=false(nt,nx);pinch=false(nt,nx);deep=false(nt,nx);deep(t>=1.05,model.live)=true;
for j=1:nx
 if any(abs(x(j)-g.faults)<=90),fault(t>=mt.cumtime(61,j),j)=true;end
 if x(j)>=g.wedge_tip-450&&x(j)<=g.wedge_tip+450
  iz=find(mt.strata(:,j)>=g.wedge_depth,1);if ~isempty(iz),pinch(abs(t-mt.cumtime(iz,j))<=.04,j)=true;end
 end
end
fault(:,~model.live)=false;pinch(:,~model.live)=false;
lowE=false(nt,nx);
for k=1:numel(model.datasets)
 T=model.datasets(k).target;L=sqrt(movmean(T.^2,21,1));rowrms=sqrt(mean(T.^2,2));sig=rowrms>.01*max(rowrms);
 vals=L(sig,model.live);lo=prctile(vals(:),20);hi=prctile(vals(:),45);
 m=false(nt,nx);m(sig,model.live)=L(sig,model.live)>=lo&L(sig,model.live)<=hi;lowE(:,:,k)=m;
end
end
