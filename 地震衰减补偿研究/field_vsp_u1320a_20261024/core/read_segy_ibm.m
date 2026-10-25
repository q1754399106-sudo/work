function [X,dt,H]=read_segy_ibm(path,first,count)
% Read IBM-float SEG-Y traces first..first+count-1 (1-based). H holds raw
% big-endian header words used here: seq (1-4), sp (17-20), delay ms (109-110),
% transit ms (235-236, LDEO VSP convention).
fid=fopen(path,'r','ieee-be');c=onCleanup(@()fclose(fid));
fseek(fid,3216,'bof');dt=fread(fid,1,'uint16')/1e6;fseek(fid,3220,'bof');ns=fread(fid,1,'uint16');
fseek(fid,3224,'bof');fmt=fread(fid,1,'uint16');assert(fmt==1,'IBM float only');
fseek(fid,0,'eof');ntr=(ftell(fid)-3600)/(240+4*ns);assert(ntr==round(ntr));
if nargin<2,first=1;count=ntr;end;last=min(ntr,first+count-1);n=last-first+1;
X=zeros(ns,n);H=struct('seq',zeros(1,n),'sp',zeros(1,n),'delay_ms',zeros(1,n),'transit_ms',zeros(1,n));
for k=1:n
 fseek(fid,3600+(first+k-2)*(240+4*ns),'bof');h=fread(fid,240,'uint8=>uint8')';
 w=@(o)double(typecast(fliplr(h(o:o+3)),'int32'));s=@(o)double(typecast(fliplr(h(o:o+1)),'int16'));
 H.seq(k)=w(1);H.sp(k)=w(17);H.delay_ms(k)=s(109);H.transit_ms(k)=s(235);
 u=fread(fid,ns,'uint32=>uint32');X(:,k)=ibm(u);
end
end
function v=ibm(u)
sgn=1-2*double(bitshift(u,-31));ex=double(bitand(bitshift(u,-24),uint32(127)))-64;
fr=double(bitand(u,uint32(16777215)))/16777216;v=sgn.*fr.*16.^ex;
end
