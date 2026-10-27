function score_field()
% Scoring of PLAN_ZH.md ("评价", "门槛与判定") for F-11 T2: T1 (NIVSP corridor CORR65_W1), T2 (AC/DEN
% synthetic), T1' (CORR65_W3), gates G2, G2b, G3, G4, G5 and the registered descriptive quantities.
% Uses no method output to set anything.
P=setup_paths();root=P.root;in=load(fullfile(root,'inputs','curtain.mat'));tw=in.tw;fp=in.fp;
vq=jsondecode(fileread(fullfile(root,'analysis','vsp_q.json')));raw=fullfile(root,'raw','vsp_velocity');
% T1 / T1' references: corridor stacks (SLB normal polarity -> x(-1)), same zero-phase filter as the outputs.
ni=read_table(fullfile(raw,'VSPNI_COMPUTED_39.ASC'),7);rs=read_table(fullfile(raw,'VSPRS_COMPUTED_42.ASC'),7);
cref=@(T,c)-interp1(T(:,1),zpfilter(T(:,c),.001),tw);
ref1=cref(ni,2);ref1w3=cref(ni,5);refrs=cref(rs,2);
nz=ni(ni(:,2)~=0,1);w1=tw>=max(min(nz)+.020,1.25)-1e-9&tw<=min(max(nz)-.020,2.95)+1e-9;
wrs=tw>=1.25-1e-9&tw<=2.0+1e-9;
% T2 reference: impedance DEN*304800/AC, MD -> TWT by the combined time-depth table, 1 ms bins,
% RC between consecutive bins, convolved with Ricker(fp) on 1 ms, sampled at 4 ms.
lg=readmatrix(fullfile(root,'inputs','logs_md.csv'));cs=read_table(fullfile(raw,'TZV_TIME_CHECKSHOT_3.ASC'),3);
[~,u]=unique(cs(:,1));cs=cs(u,:);ok=all(isfinite(lg(:,2:3)),2);md=lg(ok,1);Z=lg(ok,3)*304800./lg(ok,2);
tz=interp1(cs(:,1),cs(:,3),md)/1000;ok=isfinite(tz);tz=tz(ok);Z=Z(ok);
t1ms=(floor(min(tz)*1000):ceil(max(tz)*1000))'/1000;bin=round((tz-t1ms(1))*1000)+1;
lz=accumarray(bin,log(Z),[numel(t1ms) 1],@mean,NaN);have=isfinite(lz);lz=interp1(t1ms(have),lz(have),t1ms);
Zb=exp(lz);rc=[0;diff(Zb)./(Zb(1:end-1)+Zb(2:end))];
k=(-100:100)'/1000;z=pi*fp*k;ric=(1-2*z.^2).*exp(-z.^2);syn=conv(rc,ric,'same');
ref2=interp1(t1ms,syn,tw);w2=tw>=max(t1ms(1)+.020,1.25)-1e-9&tw<=min(t1ms(end)-.020,2.95)+1e-9;
% Scores.
methods=[{'Input'} P.methods {'SC','QSC'}];rows={};runs={P.runs.id};R=struct();
for i=1:numel(runs)
 L=load(fullfile(root,'results',[runs{i} '.mat']));R.(runs{i})=L;
 L.G.SC=sc(root,runs{i},'SC');L.G.QSC=sc(root,runs{i},'QSC');
 for m=methods
  A=L.G.(m{1});x=path_trace(A,in.idx);xf=zpfilter(x,.004);
  [t1,s1]=tie(xf,ref1,tw,w1);[t1b,s1b]=tie(xf,ref1w3,tw,w1);[t2,s2]=tie(x,ref2,tw,w2);[trs,srs]=tie(xf,refrs,tw,wrs);
  [hf,band]=spectral(A,L.G.Input,tw);
  rows(end+1,:)={runs{i},L.r.e,m{1},t1,s1,t1b,s1b,t2,s2,trs,srs,hf,band(1),band(2)}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'run','q_error_percent','method','T1','T1_shift_ms','T1w3','T1w3_shift_ms', ...
 'T2','T2_shift_ms','T1rs','T1rs_shift_ms','hf_ratio_60_100Hz_2to3s','band05_Hz','band95_Hz'});
an=fullfile(root,'analysis');writetable(T,fullfile(an,'scores.csv'));
v=@(run,m,col)T.(col)(strcmp(T.run,run)&strcmp(T.method,m));
verdict=@(d,a,b)ifelse(d>=.02,a,ifelse(d<=-.02,b,'no discernible difference'));
% G2, G2b (T1 primary, T2 secondary), G3.
d2=v('V0','SC','T1')-v('V0','QSC','T1');d2t2=v('V0','SC','T2')-v('V0','QSC','T2');
d2b=v('V0','SC','T1')-v('V0','Standard','T1');d2bt2=v('V0','SC','T2')-v('V0','Standard','T2');
in1=v('V0','Input','T1');
% G4: implied scale against the 90% interval of Q_NI.
lo=vq.NI.Q_90CI(1)/vq.NI.Q;hi=vq.NI.Q_90CI(2)/vq.NI.Q;q=struct('run',{},'start_scale',{},'selected_factor',{},'implied_scale',{},'inside',{});
for id={'Vm50','Vm30','V0','Vp30','Vp50'}
 L=R.(id{1});s=(1+L.r.e/100)*L.info.qscan_factor;q(end+1)=struct('run',id{1},'start_scale',1+L.r.e/100, ...
  'selected_factor',L.info.qscan_factor,'implied_scale',s,'inside',s>=lo&&s<=hi); %#ok<AGROW>
end
% G5: QSC - SC at the four wrong starts.
wrong={'Vm50','Vm30','Vp30','Vp50'};g5=cellfun(@(id)v(id,'QSC','T1')-v(id,'SC','T1'),wrong);
g5t2=cellfun(@(id)v(id,'QSC','T2')-v(id,'SC','T2'),wrong);
g5v=ifelse(sum(g5>=.02)>=3,'supports the rule (QScan+SC better)',ifelse(sum(g5<=-.02)>=3,'opposite (SC better)','no discernible difference'));
% Input constant-phase scan (descriptive).
xin=path_trace(R.V0.G.Input,in.idx);ph=-90:5:90;pt=zeros(size(ph));H=taper_pct('hilbert',xin);
for p=1:numel(ph),pt(p)=tie(zpfilter(cosd(ph(p))*xin-sind(ph(p))*H,.004),ref1,tw,w1);end
[pmax,pi_]=max(pt);
G=struct('G1_all_pass',vq.G1.all_pass,'G1_RS',vq.G1.RS.G1_pass,'G1_NI',vq.G1.NI.G1_pass,'G1b_pass',vq.G1b_pass, ...
 'descriptive_only',~vq.G1.all_pass, ...
 'G2_delta_T1',d2,'G2_verdict_T1',verdict(d2,'SC better','QScan+SC better'),'G2_delta_T2',d2t2, ...
 'G2_verdict_T2',verdict(d2t2,'SC better','QScan+SC better'), ...
 'G2b_delta_T1',d2b,'G2b_verdict_T1',verdict(d2b,'SC better','Standard better'),'G2b_delta_T2',d2bt2, ...
 'G2b_verdict_T2',verdict(d2bt2,'SC better','Standard better'), ...
 'G3_SC',v('V0','SC','T1')>=in1+.02,'G3_QSC',v('V0','QSC','T1')>=in1+.02,'G3_Standard',v('V0','Standard','T1')>=in1+.02, ...
 'G4_pass',sum([q.inside])>=4,'G4_interval_scale',[lo hi],'qscan',q, ...
 'G5_QSC_minus_SC_T1',g5,'G5_verdict',g5v,'G5_QSC_minus_SC_T2',g5t2, ...
 'polarity_flag',in1<0||v('V0','Input','T2')<0,'input_T1',in1,'input_T2',v('V0','Input','T2'), ...
 'input_T1_minus_T1w3',in1-v('V0','Input','T1w3'),'input_phase_scan_best_deg',ph(pi_),'input_phase_scan_best_T1',pmax, ...
 'T1_window_s',[min(tw(w1)) max(tw(w1))],'T2_window_s',[min(tw(w2)) max(tw(w2))],'Q_VSP',reshape(vq.Q_models.VSP,1,[]),'fpeak',fp);
uq_writejson(fullfile(an,'gates.json'),G);disp(T(:,[1 3 4 6 8 10]));disp(jsonencode(G));
save(fullfile(an,'references.mat'),'tw','ref1','ref1w3','refrs','ref2','w1','w2','wrs','ph','pt','t1ms','rc');
end

function X=sc(root,run,label)
L=load(fullfile(root,'results_sc',sprintf('F11_%s_%s.mat',run,label)));X=L.X;
end
function x=path_trace(A,idx)
x=zeros(size(A,1),1);for n=1:size(A,1),x(n)=mean(A(n,idx(n)-2:idx(n)+2));end
end
function [best,lag]=tie(x,ref,tw,w)
best=-Inf;lag=NaN;tf=tw(w);rr=ref(w);
for s=-.020:.001:.020
 y=interp1(tw,x,tf+s,'linear',NaN);ok=isfinite(y)&isfinite(rr);a=y(ok)-mean(y(ok));b=rr(ok)-mean(rr(ok));
 c=(a'*b)/sqrt((a'*a)*(b'*b));if c>best,best=c;lag=s*1000;end
end
end
function y=zpfilter(x,dt)
% Zero-phase: |H| = 3rd-order Butterworth high-pass 6 Hz x 12th-order low-pass 65 Hz.
x=x(:);n=numel(x);N=2^nextpow2(2*n);f=abs([0:N/2 -N/2+1:-1]')/(N*dt);
H=1./sqrt(1+(6./max(f,realmin)).^6)./sqrt(1+(f/65).^24);H(1)=0;y=real(ifft(fft(x,N).*H));y=y(1:n);
end
function [r,band]=spectral(A,I,tw)
w=tw>=2.0-1e-9&tw<=3.0+1e-9;n=nnz(w);f=(0:n-1)'/(n*.004);b=f>=60&f<=100;
FA=abs(fft(A(w,:))).^2;FI=abs(fft(I(w,:))).^2;r=sum(sum(FA(b,:)))/sum(sum(FI(b,:)));
h=f<=125;p=cumsum(sum(FA(h,:),2));p=p/p(end);fh=f(h);band=[fh(find(p>=.05,1)) fh(find(p>=.95,1))];
end
function c=read_table(path,ncol)
% Rows with exactly ncol numbers from a contractor ASCII listing.
L=splitlines(fileread(path));c=zeros(0,ncol);
for k=1:numel(L),v=sscanf(L{k},'%f')';if numel(v)==ncol,c(end+1,:)=v;end,end %#ok<AGROW>
end
function o=ifelse(c,a,b),if c,o=a;else,o=b;end,end
