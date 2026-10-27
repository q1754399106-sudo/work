function prep_line()
% Data-driven preparation shared by all methods (PLAN_ZH.md "测线准备"; experiment 5 prep_line.m with the
% receivers of segments RS and NI): 99-trace curtain along the well path (bilinear from the 25 m bins),
% window 1.200-3.000 s, dominant frequency fp, blind Haar sigma, and the per-sample scoring trace.
P=setup_paths();root=P.root;vq=jsondecode(fileread(fullfile(root,'analysis','vsp_q.json')));
geo=jsondecode(fileread(fullfile(root,'probe','vsp_geometry.json')));
rs=asrow(geo.RS);rs=rs([rs.md]<=P.rs_md_max);lv=[rs asrow(geo.NI)];
md=[lv.md]';rx=[lv.rec_x]';ry=[lv.rec_y]';fb=[lv.fb_ms]'/1000;nl=numel(lv);
% Polyline by MD (shallow -> deep), arc length s; 99 points at 25 m centred on the arc midpoint.
[~,o]=sort(md);px=rx(o);py=ry(o);s=[0;cumsum(hypot(diff(px),diff(py)))];[s,u]=unique(s,'stable');px=px(u);py=py(u);
send=s(end);sk=send/2+25*((1:99)'-50);cx=zeros(99,1);cy=cx;
at=@(q)[interp1(s,px,q) interp1(s,py,q)];
d0=at(0)-at(100);d0=d0/norm(d0);d1=at(send)-at(send-100);d1=d1/norm(d1);
for k=1:99
 if sk(k)<0,p=at(0)+(-sk(k))*d0;elseif sk(k)>send,p=at(send)+(sk(k)-send)*d1;else,p=at(sk(k));end
 cx(k)=p(1);cy(k)=p(2);
end
% Cube headers and bilinear weights on the 25 m (every 2nd) grid.
fid=fopen(P.cube,'r','ieee-be');ns=1126;tb=240+4*ns;fseek(fid,0,'eof');ntr=(ftell(fid)-3600)/tb;fseek(fid,3600,'bof');
raw=fread(fid,[tb ntr],'uint8=>uint8');fclose(fid);
il=w32(raw,189);xl=w32(raw,193);sc=w16(raw,71);sc(sc<0)=-1./sc(sc<0);sc(sc==0)=1;bx=w32(raw,181).*sc;by=w32(raw,185).*sc;
A=[ones(ntr,1) bx by]\[il xl];
key=@(i,j)(i-9985)/2*1000+(j-1932)/2;map=containers.Map(num2cell(key(il,xl)),num2cell((1:ntr)'));
F=[ones(99,1) cx cy]*A;i0=9985+2*floor((F(:,1)-9985)/2);j0=1932+2*floor((F(:,2)-1932)/2);
a=(F(:,1)-i0)/2;b=(F(:,2)-j0)/2;tr=zeros(99,4);wt=tr;
for k=1:99
 c=[i0(k) j0(k);i0(k)+2 j0(k);i0(k) j0(k)+2;i0(k)+2 j0(k)+2];
 tr(k,:)=arrayfun(@(q)map(key(c(q,1),c(q,2))),1:4);
 wt(k,:)=[(1-a(k))*(1-b(k)) a(k)*(1-b(k)) (1-a(k))*b(k) a(k)*b(k)];
end
X=zeros(ns,99);
for k=1:99,for q=1:4,X(:,k)=X(:,k)+wt(k,q)*ibm(raw(241:end,tr(k,q)));end,end
t=(0:ns-1)'*.004;win=t>=1.2-1e-9&t<=3.0+1e-9;tw=t(win);D=X(win,:);assert(numel(tw)==451);
% Dominant frequency: Ricker x dictionary absorption at 1.40 s (Q_VSP), fitted to the 1.20-1.60 s spectrum.
Qv=reshape(vq.Q_models.VSP,1,[]);cm=struct('dt',.004,'nt',numel(tw),'nfft',2048,'fref',100,'fpeak',NaN,'t',tw, ...
 'vp',reshape(vq.vp,1,[]),'h',reshape(vq.h,1,[]),'Q',Qv,'solver_offset',0,'lambda_x',.001,'lambda_fixed',.001);
a14=ray_absorption(1.40,0,cm);e=t>=1.2-1e-9&t<=1.6+1e-9;nf=512;
S=mean(abs(fft(X(e,:).*taper_pct('taper',nnz(e),.2),nf)),2);f=(0:nf/2)'/(nf*.004);S=S(1:nf/2+1);bd=f>=5&f<=80;
grid=10:.25:80;mis=zeros(size(grid));
for q=1:numel(grid),r=(f(bd).^2/grid(q)^2).*exp(-f(bd).^2/grid(q)^2).*exp(-pi*f(bd)*a14);mis(q)=norm(S(bd)-(r\S(bd))*r);end
[~,q]=min(mis);fp=grid(q);cm.fpeak=fp;
sigma=zeros(1,99);for j=1:99,ni=estimate_trace_snr(D(:,j),tw,struct('method','haar'));sigma(j)=ni.sigma;end
% Scoring trace per sample: well arc length at one-way time t/2, nearest curtain trace (clamped to 3..97).
sl=zeros(nl,1);for k=1:nl,[~,m]=min(hypot(px-rx(k),py-ry(k)));sl(k)=s(m);end
[fs,o2]=sort(fb);[fs,~,g]=unique(fs);ss=accumarray(g,sl(o2),[],@mean);wx=accumarray(g,rx(o2),[],@mean);wy=accumarray(g,ry(o2),[],@mean);
tau=min(max(tw/2,fs(1)),fs(end));sw=interp1(fs,ss,tau);well_xy=[interp1(fs,wx,tau) interp1(fs,wy,tau)];
[~,idx]=min(abs(sw-sk'),[],2);idx=min(max(idx,3),97);
meta=struct('arc_length_m',send,'curtain_xy_first',[cx(1) cy(1)],'curtain_xy_last',[cx(99) cy(99)], ...
 'fpeak',fp,'tstar_twoway_at_1p40',a14,'sigma_median',median(sigma),'score_trace_range',[min(idx) max(idx)], ...
 'Q_VSP',Qv,'vp',cm.vp,'h',cm.h);
if ~isfolder(fullfile(root,'inputs')),mkdir(fullfile(root,'inputs'));end
save(fullfile(root,'inputs','curtain.mat'),'D','tw','cm','sigma','fp','sk','cx','cy','tr','wt','idx','well_xy','S','f','mis','-v7');
uq_writejson(fullfile(root,'inputs','curtain_meta.json'),meta);disp(jsonencode(meta));
end

function v=w32(raw,b),v=double(typecast(reshape(flipud(raw(b:b+3,:)),[],1),'int32'));v=v(:);end
function v=w16(raw,b),v=double(typecast(reshape(flipud(raw(b:b+1,:)),[],1),'int16'));v=v(:);end
function v=ibm(bytes)
u=typecast(reshape(flipud(reshape(bytes,4,[])),[],1),'uint32');
sgn=1-2*double(bitshift(u,-31));ex=double(bitand(bitshift(u,-24),uint32(127)))-64;
fr=double(bitand(u,uint32(16777215)))/16777216;v=sgn.*fr.*16.^ex;
end
function s=asrow(s),if iscell(s),s=[s{:}];end,s=s(:)';end
