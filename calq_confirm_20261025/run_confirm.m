function run_confirm()
% Exploratory SharedGhost confirmation on new models C1-C3 (inputs/PLAN_ZH.md).
% Solver side: calibration gathers and observations only; truth used for scoring.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));maxNumCompThreads(4);
policy=jsondecode(fileread(fullfile(root,'inputs','candidate_policy.json')));started=tic;rows={};qrows={};
specs={'C1',2026102511,[-.20 .010],[1100 2400],[-120 100],[2600 320 820],2200,600,1180,[1400 1300];
       'C2',2026102521,[-.12 .016;.03 .032],[700 1900],[150 -110],[1300 260 700],350,1900,1220,[2500 1150];
       'C3',2026102531,[-.16 .012;.04 .025],[1600 2700],[-100 130],[500 300 900],2700,1000,1150,[1800 1280]};
for i=1:size(specs,1)
 sp=specs(i,:);
 g=struct('seed',sp{2},'calibration_seed',sp{2}+1,'vp0',1860+20*i,'q0',32+2*i,'q_gradient',.034+.002*i,'fold_amp',[30+8*i 18],'fold_phase',.7*i, ...
  'faults',sp{4},'throws',sp{5},'steps',[640 110;1050 125;1520 135],'channel',sp{6},'wedge_start',sp{7},'wedge_tip',sp{8},'wedge_depth',sp{9}, ...
  'wedge_thickness',44,'gas',sp{10},'diffractors',[sp{4}(1) .78 .015;sp{4}(2) 1.02 .012;sp{6}(1) .72 .010],'dead',[20+5*i 70+4*i],'burst_traces',[15 50 85], ...
  'groundroll_x',600*i,'linear_t0',.5+.04*i,'linear_slope',(-1)^i*.00007,'ghost',sp{3}, ...
  'observations',struct('name',{'Moderate','Strong'},'snr_db',{14,8},'seed',{sp{2}+5,sp{2}+6},'coherent_fraction',{.18,.30},'multiple_scale',{1,1.5},'static_std_ms',{1.5,3}));
 model=gen_complex_ghost(g);p=model.public;qt=model.truth.calibration_Q;[fault,pinch,lowE,deep]=complex_regions(model,g);
 save(fullfile(root,'data',[sp{1} '.mat']),'p','-v7');tr=model.datasets;save(fullfile(root,'truth',[sp{1} '_truth.mat']),'g','qt','fault','pinch','lowE','deep','-v7');
 [e1,c1]=p1_prepare_candidate(model.calibrationD,p.offset,p,policy);[e2,c2]=prepare_sharedghost(model.calibrationD,p.offset,p);
 cms={'current',uq_inverse_config(c1,e1);'SharedGhost',uq_inverse_config(c2,e2)};
 for c=1:2
  qrows{end+1}=struct('model',sp{1},'calibration',cms{c,1},'Q',cms{c,2}.Q,'true_Q',qt,'log_bias',log(cms{c,2}.Q./qt),'fpeak',cms{c,2}.fpeak); %#ok<AGROW>
  fprintf('%s %-11s Q %s ln %s f %.2f\n',sp{1},cms{c,1},mat2str(cms{c,2}.Q,4),mat2str(log(cms{c,2}.Q./qt),3),cms{c,2}.fpeak);
 end
 fprintf('%s SharedGhost ghost tau %.0f ms g %.2f | true %s | %.0fs\n',sp{1},1e3*e2.ghost_delay,e2.ghost_amplitude,mat2str(sp{3}),toc(started));
 for k=1:2
  D=tr(k).D;s=zeros(1,99);for j=1:99,if all(D(:,j)==0),s(j)=realmin;else,ni=estimate_trace_snr(D(:,j),model.t);s(j)=ni.sigma;end,end
  T=struct('target',tr(k).target,'masks',struct('fault',fault,'pinch',pinch,'low_energy',lowE(:,:,k),'deep',deep));
  for c=1:2
   cm=cms{c,2};cm.t=cm.t(:);[S0,S]=solve(D,s,cm);
   for m={'FullCalRefit',S0;'Standard',S}'
    x=score_gather(m{2},T,cm);rows{end+1}=struct('model',sp{1},'condition',tr(k).name,'calibration',cms{c,1},'method',m{1},'overall',x.overall,'deep',x.deep,'low_energy',x.low_energy,'fault',x.fault,'pinch',x.pinch); %#ok<AGROW>
    fprintf('  %s %-8s %-11s %-12s overall %.4f deep %.4f lowE %.4f\n',sp{1},tr(k).name,cms{c,1},m{1},x.overall,x.deep,x.low_energy);
   end
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis','confirm_metrics.csv'));uq_writejson(fullfile(root,'analysis','confirm_q.json'),[qrows{:}]);
res=struct();
for m={'Standard','FullCalRefit'}
 a=M(strcmp(M.calibration,'SharedGhost')&strcmp(M.method,m{1}),:);b=M(strcmp(M.calibration,'current')&strcmp(M.method,m{1}),:);
 r=a.overall./b.overall;rd=a.deep./b.deep;rl=a.low_energy./b.low_energy;
 res.(m{1})=struct('overall_ratios',r','deep_ratios',rd','low_energy_ratios',rl','mean_overall',mean(r),'mean_deep',mean(rd),'max_overall',max(r));
 fprintf('%s: overall ratio mean %.4f (max %.4f), deep mean %.4f, lowE mean %.4f\n',m{1},mean(r),max(r),mean(rd),mean(rl));
end
Q=[qrows{:}];lb=@(cal)mean(arrayfun(@(q)abs(q.log_bias(4)),Q(strcmp({Q.calibration},cal))));
s=res.Standard;gates=struct('K2',s.mean_overall<=.97,'K3',s.mean_deep<=.95,'K4',s.max_overall<=1.02,'K1_reported_layer4_reduction',1-lb('SharedGhost')/lb('current'));
fprintf('K2 %d (%.4f)  K3 %d (%.4f)  K4 %d (max %.4f)  K1 layer-4 |ln| reduction %.3f\n',gates.K2,s.mean_overall,gates.K3,s.mean_deep,gates.K4,s.max_overall,gates.K1_reported_layer4_reduction);
uq_writejson(fullfile(root,'analysis','confirm_summary.json'),struct('exploratory',true,'gates',gates,'ratios',res,'elapsed_seconds',toc(started)));
end

function [S0,S]=solve(D,s,cm)
d=prune_dictionary(cm,450);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);[nt,nx]=size(D);
R0=zeros(nt,nx);S0=R0;sup=cell(1,nx);E=sup;V=sup;
for b=1:9,ids=(b-1)*11+(1:11);ini=reflectivity_refit_prepared(D(:,ids),dicts,s(ids));so=reflectivity_soft_prepared(D(:,ids),dicts,s(ids),ini);
 R0(:,ids)=ini.reflectivity;S0(:,ids)=ini.gather;sup(ids)=ini.support;E(ids)=so.group_energy;V(ids)=so.group_noise_energy;end
cs=coherence_scores(R0,sup,d,any(D~=0,1),2,.25);st=struct('R0',R0,'supports',{sup},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'W',d.W);
S=apply_gain(st,cs.c,struct('lambda_min',.5,'lambda_max',16));
end
