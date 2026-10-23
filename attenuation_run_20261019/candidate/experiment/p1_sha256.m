function digest=p1_sha256(input)
% File bytes for char path; raw bytes for uint8. Java is bundled with MATLAB.
if ischar(input) || isstring(input)
 fid=fopen(input,'rb'); assert(fid>=0,'Cannot hash file: %s',input);
 cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
 bytes=fread(fid,Inf,'*uint8');
else, bytes=uint8(input(:)); end
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
