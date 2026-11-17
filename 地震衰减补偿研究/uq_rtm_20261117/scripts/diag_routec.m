function diag_routec()
% Diagnostic (no scores): route-C filtered traces against their input on one shot.
P=rtm_setup();D=load(fullfile(P.cache,'data_P1.mat'));k=26;
for tg={'m12','C_m100','att'}
 X=squeeze(D.(tg{1})(:,:,k));e=sum(abs(X).^2,1);
 fprintf('%-7s energy per freq (every 15th): %s\n',tg{1},sprintf('%9.2e',e(1:15:end)));
end
a=squeeze(D.m12(:,:,k));c=squeeze(D.C_m100(:,:,k));
for j=1:15:numel(P.f)
 r=(a(:,j)'*c(:,j))/norm(a(:,j))/norm(c(:,j));fprintf('f %4.1f Hz: |corr(input,C)| %.3f phase %.2f\n',P.f(j),abs(r),angle(r));
end
end
