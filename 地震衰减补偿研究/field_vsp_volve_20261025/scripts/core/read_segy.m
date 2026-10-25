function [X,dt,Hb]=read_segy(path,idx)
% IBM-float SEG-Y reader. idx: 1-based trace numbers (all if omitted). Hb: 240-by-n raw header bytes.
fid=fopen(path,'r','ieee-be');c=onCleanup(@()fclose(fid));
fseek(fid,3216,'bof');dt=fread(fid,1,'uint16')/1e6;fseek(fid,3220,'bof');ns=fread(fid,1,'uint16');
fseek(fid,3224,'bof');fmt=fread(fid,1,'uint16');assert(fmt==1,'IBM float only');
tb=240+4*ns;fseek(fid,0,'eof');ntr=(ftell(fid)-3600)/tb;assert(ntr==round(ntr));
if nargin<2,idx=1:ntr;end
X=zeros(ns,numel(idx));Hb=zeros(240,numel(idx),'uint8');
for k=1:numel(idx)
 fseek(fid,3600+(idx(k)-1)*tb,'bof');Hb(:,k)=fread(fid,240,'uint8=>uint8');
 X(:,k)=ibm(fread(fid,ns,'uint32=>uint32'));
end
end
function v=ibm(u)
sgn=1-2*double(bitshift(u,-31));ex=double(bitand(bitshift(u,-24),uint32(127)))-64;
fr=double(bitand(u,uint32(16777215)))/16777216;v=sgn.*fr.*16.^ex;
end
