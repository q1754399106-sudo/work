function make_sealed_data()
% Sealed V4/V5 generation after FROZEN.json; solver inputs and truth kept apart.
root=fileparts(mfilename('fullpath'));base=fileparts(root);v1=fullfile(base,'coherence_soft_20261023');
addpath(fullfile(root,'core'),fullfile(v1,'core'),fullfile(v1,'source'));maxNumCompThreads(4);
assert(isfile(fullfile(root,'inputs','FROZEN.json')));out=fullfile(root,'inputs','sealed');truth=fullfile(root,'truth');
assert(~isfolder(out)&&~isfolder(truth),'Sealed data exist; refusing to regenerate.');mkdir(out);mkdir(truth);
primary=jsondecode(fileread(fullfile(base,'soft_refit_20261023','inputs','primary.json')));
cm=primary.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';cm.Q=[38 70 60 140];nx=99;nt=cm.nt;dt=cm.dt;t=cm.t;u=linspace(0,10,nx);
times0=[.20 .29 .41 .55 .69 .86 1.01 1.21 1.47];amps=[.08 .009 -.06 .010 .07 -.008 -.045 .006 .035];
R=zeros(nt,nx);ev=nan(nx,9);fault=false(nt,nx);pinch=false(nt,nx);
for j=1:nx
 ti=times0-.025*(j-1)/98+.005*sin(u(j)*.7+(0:8)*1.1);
 if j>=35,ti(times0>=.50)=ti(times0>=.50)-.020;end
 if j>=75,ti(times0>=.80)=ti(times0>=.80)+.016;end
 sep=interp1([1 25 60 99],[.03 .03 .14 .14],j);ti(4)=ti(3)+sep;
 a=amps.*(1+.10*sin(u(j)*.5+(0:8)*.7));a(4)=a(4)*min(max((j-25)/20,0),1);
 ii=round(ti/dt)+1;ev(j,:)=ii;ev(j,a==0)=NaN;
 for e=1:9,if a(e)~=0,R(ii(e),j)=R(ii(e),j)+a(e);end,end
 if (j>=31&&j<=39),fault(t>=.50,j)=true;end
 if (j>=71&&j<=79),fault(t>=.80,j)=true;end
 if j<=45,pinch(abs(t-ti(4))<=.04,j)=true;end
end
cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{2026102441,2026102442});ds=gen_layered_scene(cm,R,cond);
config=cm;obs=rmfield(ds,'target');save(fullfile(out,'V4.mat'),'config','obs','-v7');
tr=struct('name',{ds.name},'target',{ds.target});regions=struct('fault',fault,'pinch',pinch);reflectivity=R;events=ev;weak=[2 4 6 8]; %#ok<NASGU>
save(fullfile(truth,'V4_truth.mat'),'tr','reflectivity','events','weak','regions','-v7');
g=struct('seed',2026102451,'calibration_seed',2026102452,'vp0',1920,'q0',36,'q_gradient',.032, ...
 'fold_amp',[52 15],'fold_phase',2.2,'faults',[1500 2550],'throws',[140 -95],'steps',[700 120;1150 110;1600 150], ...
 'channel',[700 280 880],'wedge_start',2800,'wedge_tip',1300,'wedge_depth',1250,'wedge_thickness',40,'gas',[2300 1100], ...
 'diffractors',[1500 .82 .015;2550 1.05 .012;700 .78 .010],'dead',[18 84],'burst_traces',[25 60 93], ...
 'groundroll_x',300,'linear_t0',.62,'linear_slope',.00006, ...
 'observations',struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{2026102461,2026102462}, ...
  'coherent_fraction',{.18,.30},'multiple_scale',{1,1.5},'static_std_ms',{1.5,3}));
started=tic;model=gen_complex_model(g);policy=jsondecode(fileread(fullfile(v1,'inputs','candidate_policy.json')));p=model.public;
[est,corrected,cal]=p1_prepare_candidate(model.calibrationD,p.offset,p,policy);config=uq_inverse_config(corrected,est);
fprintf('V5 calibrated %.1fs: Q %s f %.2f, coarse true %s\n',toc(started),mat2str(config.Q,4),config.fpeak,mat2str(model.truth.calibration_Q,4));
[fault,pinch,lowE,deep]=complex_regions(model,g);obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
for k=1:2
 D=model.datasets(k).D;sigma=zeros(1,nx);
 for j=1:nx,if all(D(:,j)==0),sigma(j)=realmin;else,ni=estimate_trace_snr(D(:,j),model.t);sigma(j)=ni.sigma;end,end
 obs(k)=struct('name',model.datasets(k).name,'D',D,'sigma',sigma);tr(k)=struct('name',model.datasets(k).name,'target',model.datasets(k).target);
end
calibration=struct('estimate',est,'diagnostics',cal); %#ok<NASGU>
save(fullfile(out,'V5.mat'),'config','obs','calibration','-v7');
regions=struct('fault',fault,'pinch',pinch,'low_energy',lowE,'deep',deep);truthQ=model.truth.calibration_Q; %#ok<NASGU>
save(fullfile(truth,'V5_truth.mat'),'tr','regions','g','truthQ','-v7');
end
