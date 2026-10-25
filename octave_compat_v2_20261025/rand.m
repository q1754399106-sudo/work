function X=rand(varargin)
% Octave shim (v2): MATLAB rand after rng(seed,'twister'): genrand_res53 of two words per value.
global ML2
dims=ml2_dims(varargin);n=prod(dims);ml2_words(2*n);k=ML2.k;
a=ML2.W(k:2:k+2*n-2);b=ML2.W(k+1:2:k+2*n-1);ML2.k=k+2*n;
X=reshape((floor(a/32)*67108864+floor(b/64))/9007199254740992,dims);
end
