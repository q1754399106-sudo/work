function X=randn(varargin)
% Octave shim (v2): MATLAB randn after rng(seed,'twister') (256-strip Ziggurat, see mlrandn.py),
% implemented in Octave on the buffered word stream shared with rand.m.
global ML2
persistent T F R
if isempty(T)
 f=@(x) exp(-.5*x*x);dn=3.654152885361009;vn=0.004928673233974658;T=zeros(1,257);T(1)=vn/f(dn);T(2)=dn;
 for i=3:256,T(i)=sqrt(-2*log(vn/T(i-1)+f(T(i-1))));end
 F=exp(-.5*T.^2);R=T(2:257)./T(1:256);   % index i (python) -> i+1
end
dims=ml2_dims(varargin);n=prod(dims);X=zeros(n,1);m=0;
while m<n
 ml2_words(64+3*(n-m));W=ML2.W;k=ML2.k;lim=numel(W)-16;
 while m<n&&k<lim
  a=W(k);b=W(k+1);k=k+2;i=255-floor(b/16777216);
  u=2*((floor(a/8)*16777216+mod(b,16777216))/9007199254740992)-1;x=u*T(i+1);
  if i==0
   if abs(u)<R(1),m=m+1;X(m)=x;continue;end
   while true
    xx=-log((floor(W(k)/32)*67108864+floor(W(k+1)/64))/9007199254740992)/T(2);
    yy=-log((floor(W(k+2)/32)*67108864+floor(W(k+3)/64))/9007199254740992);k=k+4;
    if yy+yy>=xx*xx,break;end
   end
   m=m+1;if u>0,X(m)=T(2)+xx;else,X(m)=-(T(2)+xx);end;continue
  end
  if abs(u)<R(i+1),m=m+1;X(m)=x;continue;end
  U=(floor(W(k)/32)*67108864+floor(W(k+1)/64))/9007199254740992;k=k+2;
  if F(i+1)+U*(F(i+2)-F(i+1))<exp(-.5*x*x),m=m+1;X(m)=x;end
 end
 ML2.k=k;
end
X=reshape(X,dims);
end
