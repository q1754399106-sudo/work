function s=gain_s(g_db)
% s of the stabilised inverse-Q amplitude (b+s)/(b^2+s) whose maximum gain is g_db (semi.gain_sigma2).
target=10^(g_db/20);lo=-40;hi=0;
for it=1:200
 s=exp((lo+hi)/2);b=-s+sqrt(s*s+s);g=(b+s)/(b*b+s);
 if g>target,lo=(lo+hi)/2;else,hi=(lo+hi)/2;end
end
s=exp((lo+hi)/2);
end
