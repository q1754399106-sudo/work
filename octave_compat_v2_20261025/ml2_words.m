function ml2_words(need)
% Make sure the buffer holds at least ML2.k-1+need words (reloads from the stream start).
global ML2
assert(~isempty(ML2),'call rng(seed,''twister'') first.');
if numel(ML2.W)>=ML2.k-1+need,return;end
n=max(2*numel(ML2.W),ML2.k-1+need+2^20);here=fileparts(mfilename('fullpath'));out=[tempname() '.bin'];
[s,t]=system(sprintf('python3 -I "%s" %d %d "%s"',fullfile(here,'mlwords.py'),ML2.seed,n,out));assert(s==0,t);
fid=fopen(out,'r');ML2.W=fread(fid,n,'uint32=>double',0,'ieee-le');fclose(fid);delete(out);
end
