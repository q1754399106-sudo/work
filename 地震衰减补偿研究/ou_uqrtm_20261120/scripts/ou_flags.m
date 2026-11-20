function F=ou_flags(M,P)
% Interpreter mask (full grid, logical) where the local prior std is s_flag instead of s_bg
% (PLAN_ZH.md section 3). Built only from the geometry an interpreter would draw on the
% migration velocity / structure: lenses and chimneys (outline), a 40 m band around a fault
% plane, and layers thinner than 100 m. No Q value of the true model is used.
[X,Z]=meshgrid(P.x,P.z);Zc=max(Z,0);S=M.spec;F=false(size(X));
if isfield(S,'lens')
 for i=1:size(S.lens,1)
  L=S.lens(i,:);F=F|((X-L(1))/L(3)).^2+((Zc-L(2))/L(4)).^2<=1;
 end
end
if isfield(S,'chimney')&&~isempty(S.chimney)
 C=S.chimney;F=F|(abs(X-C(1))<=C(2)&Zc>=C(3)&Zc<=C(4));
end
if isfield(S,'fault')&&~isempty(S.fault)
 f=S.fault;d=abs((X-f(1))*sin(f(2)*pi/180)-Zc*cos(f(2)*pi/180));F=F|(d<=40&Zc>S.z_top);
end
tops=[{@(x) S.z_top+0*x} S.iface];
for k=1:numel(tops)-1
 th=tops{k+1}(P.x)-tops{k}(P.x);
 if max(th)<100,F=F|(Zc>=tops{k}(X)&Zc<tops{k+1}(X));end
end
end
