function well_tie(line)
% E4-2/E4-3 scoring: QScan Q recovery and log-synthetic tie of every method output.
root=fileparts(mfilename('fullpath'));bae=fullfile(fileparts(root),'bae_q_20261024');
addpath(fullfile(root,'core'),fullfile(bae,'source'));
vq=jsondecode(fileread(fullfile(root,'analysis','vsp_q.json')));in=load(fullfile(root,'inputs',[line '.mat']));cm=in.cm;tw=cm.t;
std=fullfile(root,'raw','ldeo','standard');[zv,vp]=read_ldeo_log(fullfile(std,'308-U1320A_dsi2.dat'),'VCO');[zr,rho]=read_ldeo_log(fullfile(std,'308-U1320A_hldsm.dat'),'RHOM');
z=(max(min(zv),min(zr)):.1524:min(max(zv),max(zr)))';v=fill(interp1(zv,vp,z));r=fill(interp1(zr,rho,z));
v=movmedian(v,7);r=movmedian(r,7);Z=v.*r;
% Time-depth: checkshot one-way transit; above 87 mbsf integrate log velocity, 0-62 mbsf at 1.55 km/s.
dep=vq.depth_mbsf(:);tau=vq.transit_s(:);[dep,u]=unique(dep);tau=tau(u);
above=z<dep(1);t87=tau(1);tz=zeros(size(z));tz(~above)=interp1(dep,tau,z(~above),'linear','extrap');
dz=.1524;cum=cumsum(dz./(v(above)*1000),'reverse');tz(above)=t87-cum;
tau_sf=tz(1)-z(1)/1550;twt=in.meta.seafloor_twt+2*(tz-tau_sf);
tf=(twt(1):.0005:twt(end))';Zf=interp1(twt,Z,tf);rf=[0;diff(Zf)./(Zf(2:end)+Zf(1:end-1))];
fp=cm.fpeak;k=(-80:80)'*.0005;q=pi*fp*k;ric=(1-2*q.^2).*exp(-q.^2);syn=conv(rf,ric,'same');
files=dir(fullfile(root,'results',[line '_Q*.mat']));rows={};
for f=reshape(files,1,[])
 e=sscanf(regexprep(f.name,'^.*_Q([+-]\d+)\.mat$','$1'),'%d');L=load(fullfile(f.folder,f.name));
 for m=reshape(fieldnames(L.g),1,[])
  A=L.g.(m{1});tr=mean(A(:,48:52),2);best=-Inf;lag=0;
  for s=-.020:.001:.020
   x=interp1(tw,tr,tf+s,'linear',NaN);ok=isfinite(x);c=corr_(x(ok),syn(ok));if c>best,best=c;lag=s;end
  end
  hi=hf_ratio(A,L.g.Input,tw,cm.dt);
  rows(end+1,:)={line,e,m{1},best,lag*1000,hi}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'line','Q_error_percent','method','tie_corr','tie_shift_ms','hf_energy_ratio_60_100Hz'});
T=sortrows(T,{'method','Q_error_percent'});writetable(T,fullfile(root,'analysis',[line '_tie.csv']));
q=struct('start_factor',{},'selected_factor',{},'implied_Q',{},'in_vsp_ci',{});
for e=[-50 -30 0 40 60]
 L=load(fullfile(root,'results',sprintf('%s_Q%+d.mat',line,e)),'info');iq=L.info.qscan.applied_Q(2);
 q(end+1)=struct('start_factor',1+e/100,'selected_factor',L.info.qscan.selected_factor,'implied_Q',iq,'in_vsp_ci',iq>=vq.Q_SR_90CI(1)&&iq<=vq.Q_SR_90CI(2)); %#ok<AGROW>
end
get=@(m,e)T.tie_corr(strcmp(T.method,m)&T.Q_error_percent==e);rob=@(m)min(T.tie_corr(strcmp(T.method,m)&ismember(T.Q_error_percent,[-50 -30 0 30 50])));
G=struct('G2_qscan_recovery',sum([q.in_vsp_ci])>=4,'qscan',q,'G3_full',get('FullCalRefit',0)>=get('Input',0)+.02, ...
 'G3_qscan',get('QScanRefit',0)>=get('Input',0)+.02,'input_tie',get('Input',0),'full_tie',get('FullCalRefit',0),'qscan_tie',get('QScanRefit',0), ...
 'G4_robust',rob('QScanRefit')>=rob('FullCalRefit'),'worst_tie_full',rob('FullCalRefit'),'worst_tie_qscan',rob('QScanRefit'),'vsp_gate_G1',vq.gate.G1_pass, ...
 'log_window_twt',[twt(1) twt(end)],'synthetic_fpeak',fp);
uq_writejson(fullfile(root,'analysis',[line '_gates.json']),G);disp(T);disp(jsonencode(G));
save(fullfile(root,'analysis',[line '_synthetic.mat']),'tf','syn','twt','z','Z');
end
function x=fill(x),i=isfinite(x);x(~i)=interp1(find(i),x(i),find(~i),'nearest','extrap');end
function c=corr_(a,b),a=a-mean(a);b=b-mean(b);c=(a'*b)/sqrt((a'*a)*(b'*b));end
function r=hf_ratio(A,I,t,dt)
w=t>=2.6&t<=3.6;n=nnz(w);f=(0:n-1)'/(n*dt);b=f>=60&f<=100;
FA=abs(fft(A(w,:))).^2;FI=abs(fft(I(w,:))).^2;r=sum(sum(FA(b,:)))/sum(sum(FI(b,:)));
end
