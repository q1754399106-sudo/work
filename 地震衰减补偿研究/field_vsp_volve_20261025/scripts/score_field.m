function score_field()
% Scoring of PLAN_ZH.md ("评价", "门槛与判定"): T1 corridor-stack tie, T2 log-synthetic tie,
% gates G2-G4 and the registered descriptive quantities. Uses no method output to set anything.
P=setup_paths();root=P.root;in=load(fullfile(root,'inputs','curtain.mat'));tw=in.tw;fp=in.fp;
vq=jsondecode(fileread(fullfile(root,'analysis','vsp_q.json')));
cs=read_las(fullfile(root,'raw','vsp_computed','VSPNI_COMPUTED_VSPCS.LAS'));
rc=read_las(fullfile(root,'raw','vsp_velocity','TZV_TIME_SYNSEIS_1.LAS'));
% References on the 4 ms grid: corridor stack filtered like the outputs and sign-flipped (NPOL); RC * Ricker(fp).
c70=fill0(cs.VSPCS_70);ms=cs.TIME;c70f=zpfilter(c70,.001);ref1=-interp1(ms/1000,c70f,tw);
r=fill0(rc.RC);k=(-100:100)'/1000;z=pi*fp*k;ric=(1-2*z.^2).*exp(-z.^2);syn=conv(r,ric,'same');ref2=interp1(rc.TIME/1000,syn,tw);
w1=tw>=1.25-1e-9&tw<=2.65+1e-9;w2=tw>=2.29-1e-9&tw<=2.68+1e-9;
rows={};runs={P.runs.id};R=struct();
for i=1:numel(runs)
 L=load(fullfile(root,'results',[runs{i} '.mat']));R.(runs{i})=L;
 for m=[{'Input'} P.methods]
  A=L.G.(m{1});x=path_trace(A,in.idx);
  [t1,s1]=tie(zpfilter(x,.004),ref1,tw,w1);[t2,s2]=tie(x,ref2,tw,w2);[hf,band]=spectral(A,L.G.Input,tw);
  rows(end+1,:)={runs{i},L.r.model,L.r.e,m{1},t1,s1,t2,s2,hf,band(1),band(2)}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'run','q_model','q_error_percent','method','T1','T1_shift_ms','T2','T2_shift_ms', ...
 'hf_ratio_60_100Hz_2to3s','band05_Hz','band95_Hz'});
an=fullfile(root,'analysis');writetable(T,fullfile(an,'scores.csv'));
v=@(run,m,col)T.(col)(strcmp(T.run,run)&strcmp(T.method,m));
% Gates.
d1=v('V0','Standard','T1')-v('V0','QScan_Standard','T1');d2=v('V0','Standard','T2')-v('V0','QScan_Standard','T2');
verdict=@(d)ifelse(d>=.02,'independent Q better',ifelse(d<=-.02,'QScan better','no discernible difference'));
lo=vq.Q_int_90CI(1)/vq.Q_int;hi=vq.Q_int_90CI(2)/vq.Q_int;q=struct('run',{},'start_scale',{},'selected_factor',{},'implied_scale',{},'inside',{});
for id={'Vm50','Vm30','V0','Vp30','Vp50'}
 L=R.(id{1});s=(1+L.r.e/100)*L.info.qscan_factor;q(end+1)=struct('run',id{1},'start_scale',1+L.r.e/100, ...
  'selected_factor',L.info.qscan_factor,'implied_scale',s,'inside',s>=lo&&s<=hi); %#ok<AGROW>
end
rescue=struct();for id={'Vm50','Vm30','Vp30','Vp50'},rescue.(id{1})=v(id{1},'QScan_Standard','T1')-v(id{1},'Standard','T1');end
wrong=cellfun(@(id)v(id,'QScan_Standard','T1'),{'Vm50','Vm30','Vp30','Vp50'});
% Input constant-phase scan (descriptive) and the R4 curtain from the VSP files (final, inverse-Q processed).
xin=path_trace(R.V0.G.Input,in.idx);ph=-90:5:90;pt=zeros(size(ph));H=taper_pct('hilbert',xin);
for p=1:numel(ph),pt(p)=tie(zpfilter(cosd(ph(p))*xin-sind(ph(p))*H,.004),ref1,tw,w1);end
[pmax,pi_]=max(pt);r4=r4_tie(root,in,ref1,tw,w1);
G=struct('G1_pass',vq.gate.G1_pass,'G1b_pass',vq.G1b_pass,'delta1_T1',d1,'G2_verdict_T1',verdict(d1),'delta2_T2',d2, ...
 'G2_verdict_T2',verdict(d2),'G3_standard',v('V0','Standard','T1')>=v('V0','Input','T1')+.02, ...
 'G3_qscan_standard',v('V0','QScan_Standard','T1')>=v('V0','Input','T1')+.02,'G4_pass',sum([q.inside])>=4, ...
 'G4_interval_scale',[lo hi],'qscan',q,'rescue_T1',rescue,'cost_without_independent_Q_T1',v('V0','Standard','T1')-min(wrong), ...
 'polarity_flag',v('V0','Input','T1')<0||v('V0','Input','T2')<0,'input_T1',v('V0','Input','T1'),'input_T2',v('V0','Input','T2'), ...
 'input_phase_scan_best_deg',ph(pi_),'input_phase_scan_best_T1',pmax,'R4_curtain_T1',r4, ...
 'READ_standard_T1',v('R0','Standard','T1'),'READ_qscan_standard_T1',v('R0','QScan_Standard','T1'), ...
 'READ_standard_T2',v('R0','Standard','T2'),'READ_qscan_standard_T2',v('R0','QScan_Standard','T2'), ...
 'Q_VSP',reshape(vq.Q_models.VSP,1,[]),'Q_READ',reshape(vq.Q_models.READ,1,[]),'fpeak',fp);
uq_writejson(fullfile(an,'gates.json'),G);disp(T(:,1:8));disp(jsonencode(G));
save(fullfile(an,'references.mat'),'tw','ref1','ref2','w1','w2','ph','pt');
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
% Zero-phase 6(18)-70(72) Hz: |H| = 3rd-order Butterworth high-pass 6 Hz x 12th-order low-pass 70 Hz.
x=x(:);n=numel(x);N=2^nextpow2(2*n);f=abs([0:N/2 -N/2+1:-1]')/(N*dt);
H=1./sqrt(1+(6./max(f,realmin)).^6)./sqrt(1+(f/70).^24);H(1)=0;y=real(ifft(fft(x,N).*H));y=y(1:n);
end
function [r,band]=spectral(A,I,tw)
w=tw>=2.0-1e-9&tw<=3.0+1e-9;n=nnz(w);f=(0:n-1)'/(n*.004);b=f>=60&f<=100;
FA=abs(fft(A(w,:))).^2;FI=abs(fft(I(w,:))).^2;r=sum(sum(FA(b,:)))/sum(sum(FI(b,:)));
h=f<=125;p=cumsum(sum(FA(h,:),2));p=p/p(end);fh=f(h);band=[fh(find(p>=.05,1)) fh(find(p>=.95,1))];
end
function t1=r4_tie(root,in,ref1,tw,w1)
[X,dt,Hb]=read_segy(fullfile(root,'raw','vsp_computed','VSPNI_COMPUTED_17.SEGY'));
sx=hword(Hb,73,'int32')/10;sy=hword(Hb,77,'int32')/10;t=(0:size(X,1)-1)'*dt;x=zeros(numel(tw),1);
for n=1:numel(tw)
 [~,j]=min(hypot(sx-in.well_xy(n,1),sy-in.well_xy(n,2)));j=min(max(j,3),size(X,2)-2);
 x(n)=mean(interp1(t,X(:,j-2:j+2),tw(n)));
end
t1=tie(zpfilter(x,.004),ref1,tw,w1);
end
function S=read_las(path)
L=splitlines(fileread(path));sec='';names={};i=0;
for k=1:numel(L)
 l=L{k};if startsWith(l,'~'),sec=l(2);if sec=='A',i=k;break;end;continue;end
 if sec=='C'&&~isempty(strtrim(l))&&~startsWith(l,'#'),names{end+1}=matlab.lang.makeValidName(strtrim(extractBefore(l,'.'))); end %#ok<AGROW>
end
A=str2num(strjoin(L(i+1:end),newline)); %#ok<ST2NM>
A(A==-999.25)=NaN;for c=1:numel(names),S.(names{c})=A(:,c);end
end
function x=fill0(x),x(~isfinite(x))=0;end
function o=ifelse(c,a,b),if c,o=a;else,o=b;end,end
