function make_scenes()
% New V-family scenes (PLAN_ZH.md "数据"). W1-W6: layered template class, built like V4 in
% coherence_soft_v2_20261024/make_sealed_data.m with randomized event times, amplitudes, dip,
% faults and pinch-out. X1-X3: complex class, gen_complex_model + the V5 calibration chain with a
% randomized geometry g. Generators come from the round3 snapshot in cache/runner (identical to
% the code used for V4/V5). Solver inputs go to models/, truth to truth/; refuses to overwrite.
root=fileparts(fileparts(mfilename('fullpath')));runner=fullfile(root,'cache','runner');
addpath(fullfile(runner,'core'),fullfile(runner,'source'));maxNumCompThreads(4);
mdir=fullfile(root,'models');tdir=fullfile(root,'truth');
assert(~isfolder(mdir)&&~isfolder(tdir),'scenes exist; refusing to regenerate.');mkdir(mdir);mkdir(tdir);
primary=jsondecode(fileread(fullfile(root,'cache','inputs','primary.json')));
policy=jsondecode(fileread(fullfile(root,'cache','inputs','candidate_policy.json')));
for k=1:6,layered(sprintf('W%d',k),2026102711+10*(k-1),primary,mdir,tdir);end
for k=1:3,complex_scene(sprintf('X%d',k),2026102781+10*(k-1),policy,mdir,tdir);end
end

function layered(name,seed,primary,mdir,tdir)
rng(seed,'twister');
cm=primary.config;cm.t=cm.t(:);cm.vp=cm.vp(:)';cm.h=cm.h(:)';
cm.Q=[35+25*rand 40+35*rand 55+30*rand 110+35*rand];nx=99;nt=cm.nt;dt=cm.dt;u=linspace(0,10,nx);
times0=.18+.06*rand+[0 cumsum(.09+.07*rand(1,8))];
amps=zeros(1,9);amps(1:2:9)=(.03+.06*rand(1,5)).*sgn(5);amps(2:2:8)=(.006+.006*rand(1,4)).*sgn(4);
dip=(.020+.015*rand)*sgn(1);w=.5+.2*rand;phi=.7+.4*rand;
nf=1+(rand<.5);jf=sort(30+randi(41,1,nf));tf=times0(2+randi(5,1,nf));throw=(.016+.014*rand(1,nf)).*sgn(nf);
p=4+2*(rand<.5);s0=.14+.02*rand;j1=20+randi(30);j2=j1+25+randi(20);appear=rand<.5;ja=j1+randi(20);
R=zeros(nt,nx);ev=nan(nx,9);
for j=1:nx
 ti=times0+dip*(j-1)/98+.005*sin(u(j)*w+(0:8)*phi);
 for f=1:nf,if j>=jf(f),ti(times0>=tf(f))=ti(times0>=tf(f))+throw(f);end,end
 ti(p)=ti(p-1)+interp1([1 j1 j2 99],[s0 s0 .03 .03],j);
 a=amps.*(1+.10*sin(u(j)*.5+(0:8)*.7));
 if appear,a(p)=a(p)*min(max((j-ja)/20,0),1);else,a(p)=a(p)*min(max((ja+20-j)/20,0),1);end
 ii=round(ti/dt)+1;ev(j,:)=ii;ev(j,a==0)=NaN;
 for e=1:9,if a(e)~=0,R(ii(e),j)=R(ii(e),j)+a(e);end,end
end
cond=struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{seed+1,seed+2});ds=gen_layered_scene(cm,R,cond);
config=cm;obs=rmfield(ds,'target');save(fullfile(mdir,[name '.mat']),'config','obs','-v7');
tr=struct('name',{ds.name},'target',{ds.target});reflectivity=R;events=ev; %#ok<NASGU>
scene=struct('times0',times0,'amps',amps,'dip',dip,'faults',jf,'fault_time',tf,'throw',throw,'pinch_event',p,'pinch_traces',[j1 j2],'appear',appear,'taper_start',ja); %#ok<NASGU>
save(fullfile(tdir,[name '_truth.mat']),'tr','reflectivity','events','scene','-v7');
fprintf('%s: Q %s, times %s, faults %s\n',name,mat2str(cm.Q,4),mat2str(times0,3),mat2str(jf));
end

function complex_scene(name,seed,policy,mdir,tdir)
rng(seed,'twister');r=@()rand;nx=99;
f1=600+800*r();f2=min(f1+700+700*r(),2750);
if r()<.5,ws=300+200*r();wt=ws+1200+300*r();else,ws=2600+200*r();wt=ws-1200-300*r();end
chx=400+2200*r();
dead=randperm(80,2)+9;burst=randperm(88,3)+5;
g=struct('seed',seed,'calibration_seed',seed+1,'vp0',1850+100*r(),'q0',30+10*r(),'q_gradient',.030+.010*r(), ...
 'fold_amp',[30+25*r() 12+12*r()],'fold_phase',2*pi*r(),'faults',[f1 f2],'throws',(90+50*rand(1,2)).*sgn(2), ...
 'steps',[600+150*r() 100+50*r();1000+200*r() 100+50*r();1450+200*r() 110+50*r()], ...
 'channel',[chx 260+60*r() 700+200*r()],'wedge_start',ws,'wedge_tip',wt,'wedge_depth',1100+200*r(),'wedge_thickness',35+15*r(), ...
 'gas',[400+2200*r() 1050+250*r()],'diffractors',[f1 .72+.12*r() .010+.006*r();f2 .95+.12*r() .010+.006*r();chx .70+.10*r() .010], ...
 'dead',sort(dead),'burst_traces',sort(burst),'groundroll_x',300+2100*(r()<.5),'linear_t0',.55+.10*r(),'linear_slope',(.00006+.00002*r())*sgn(1), ...
 'observations',struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{seed+11,seed+12}, ...
  'coherent_fraction',{.18,.30},'multiple_scale',{1,1.5},'static_std_ms',{1.5,3}));
started=tic;model=gen_complex_model(g);p=model.public;
[est,corrected,cal]=p1_prepare_candidate(model.calibrationD,p.offset,p,policy);config=uq_inverse_config(corrected,est);
fprintf('%s calibrated %.1fs: Q %s f %.2f, coarse true %s\n',name,toc(started),mat2str(config.Q,4),config.fpeak,mat2str(model.truth.calibration_Q,4));
obs=struct('name',{},'D',{},'sigma',{});tr=struct('name',{},'target',{});
for k=1:2
 D=model.datasets(k).D;sigma=zeros(1,nx);
 for j=1:nx,if all(D(:,j)==0),sigma(j)=realmin;else,ni=estimate_trace_snr(D(:,j),model.t);sigma(j)=ni.sigma;end,end
 obs(k)=struct('name',model.datasets(k).name,'D',D,'sigma',sigma);tr(k)=struct('name',model.datasets(k).name,'target',model.datasets(k).target);
end
calibration=struct('estimate',est,'diagnostics',cal); %#ok<NASGU>
save(fullfile(mdir,[name '.mat']),'config','obs','calibration','-v7');truthQ=model.truth.calibration_Q; %#ok<NASGU>
save(fullfile(tdir,[name '_truth.mat']),'tr','g','truthQ','-v7');
end

function s=sgn(n)
s=1-2*(rand(1,n)<.5);
end
