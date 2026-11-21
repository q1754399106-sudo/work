function J=va_judge(tag)
% Judgement of P1-P3 and the stop rule (PLAN_ZH.md section 6) from analysis/va_step1_D2.json, _D4.json.
% Window W_deep only. Writes analysis/va_step1_judge<tag>.json and prints the verdicts.
if nargin<1,tag='';end
here=fileparts(mfilename('fullpath'));root=fileparts(here);J=struct();ok=struct('P1',true,'P2',true,'P3',true);
for n={'D2','D4'}
 R=jsondecode(fileread(fullfile(root,'analysis',['va_step1_' n{1} tag '.json'])));T=R.deep;m=R.m(:)';
 d=T.R0-T.R0_old;p1=d>=0.20;                                              % P1
 low=m<=1;p2a=all(T.PH100(low)<T.R0);p2b=all(T.UQW1_100(low)<T.R0);p2=p2a&&p2b;   % P2
 p3=struct();
 for q={'UQW1_25','VA'}                                                  % P3
  s=T.(q{1});a=all(s>=T.R0-0.02);b=mean(s)-mean(T.UQW1_100);
  p3.(q{1})=struct('min_minus_R0',min(s)-T.R0,'all_ge_R0_minus_002',a,'mean_minus_UQW1_100',b,'mean_ge_010',b>=0.10,'pass',a&&b>=0.10);
 end
 p3ok=p3.UQW1_25.pass&&p3.VA.pass;
 J.(n{1})=struct('R0',T.R0,'R0_old',T.R0_old,'P1_diff',d,'P1',p1, ...
  'PH100_minus_R0',T.PH100-T.R0,'UQW1_100_minus_R0',T.UQW1_100-T.R0,'P2_PH100',p2a,'P2_UQW1_100',p2b,'P2',p2,'P3',p3,'P3_pass',p3ok);
 ok.P1=ok.P1&&p1;ok.P2=ok.P2&&p2;ok.P3=ok.P3&&p3ok;
 fprintf('%s: R0 %.3f (old %.3f, diff %+.3f) P1 %d | P2 PH100 %d UQW1_100 %d | P3 UQW1_25 %d VA %d\n',n{1},T.R0,T.R0_old,d,p1,p2a,p2b,p3.UQW1_25.pass,p3.VA.pass);
end
J.P1=ok.P1;J.P2=ok.P2;J.P3=ok.P3;J.diagnosis_void=~(ok.P1&&ok.P2);
fprintf('P1 %d  P2 %d  P3 %d  -> diagnosis void (stop rule): %d\n',J.P1,J.P2,J.P3,J.diagnosis_void);
fid=fopen(fullfile(root,'analysis',['va_step1_judge' tag '.json']),'w');fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true));fclose(fid);
end
