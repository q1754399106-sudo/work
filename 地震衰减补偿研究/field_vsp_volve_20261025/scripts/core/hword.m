function v=hword(Hb,byte,type)
% Big-endian trace-header word starting at 1-based byte (type 'int32' or 'int16') for every column of Hb.
n=4;if strcmp(type,'int16'),n=2;end
v=zeros(1,size(Hb,2));
for k=1:size(Hb,2),v(k)=double(typecast(flipud(Hb(byte:byte+n-1,k)),type));end
end
