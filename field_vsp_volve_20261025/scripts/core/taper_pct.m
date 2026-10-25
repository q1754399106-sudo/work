function varargout=taper_pct(mode,varargin)
% Toolbox-free helpers. taper_pct('taper',n,r): Tukey window (cosine fraction r).
% taper_pct('pct',x,p): percentiles by linear interpolation of the sorted sample.
% taper_pct('hilbert',x): imaginary part of the analytic signal (column-wise).
switch mode
 case 'taper'
  n=varargin{1};r=varargin{2};w=ones(n,1);m=floor(r*(n-1)/2);
  if m>0,k=(0:m)';w(1:m+1)=.5*(1-cos(pi*k/m));w(end-m:end)=flipud(w(1:m+1));end
  varargout{1}=w;
 case 'pct'
  x=sort(varargin{1}(:));p=varargin{2};q=1+(numel(x)-1)*p(:)'/100;lo=floor(q);hi=ceil(q);
  varargout{1}=x(lo)'+(q-lo).*(x(hi)'-x(lo)');
 case 'hilbert'
  x=varargin{1};n=size(x,1);X=fft(x);h=zeros(n,1);h(1)=1;
  if mod(n,2)==0,h(2:n/2)=2;h(n/2+1)=1;else,h(2:(n+1)/2)=2;end
  varargout{1}=imag(ifft(X.*h));
end
end
