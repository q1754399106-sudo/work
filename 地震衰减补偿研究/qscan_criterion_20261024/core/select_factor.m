function [k,score]=select_factor(st,criterion)
% Registered candidate criteria (PLAN_ZH.md), all summed over live traces.
% A0 original QScan lasso objective; A1 lasso residual only; A2 BIC on the
% debiased refit (RSS/sigma^2 + k log n); A3 AIC on the refit (RSS/sigma^2 + 2k).
L=st.live;
switch criterion
 case 'A0',score=sum(st.J_lasso(:,L),2);
 case 'A1',score=sum(st.RSS_lasso(:,L),2);
 case 'A2',score=sum(st.RSS_refit(:,L)+st.support(:,L)*log(st.n),2);
 case 'A3',score=sum(st.RSS_refit(:,L)+2*st.support(:,L),2);
 otherwise,error('unknown criterion %s',criterion);
end
[~,k]=min(score);
end
