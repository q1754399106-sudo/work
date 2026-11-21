function J=fe2_judge(tag)
% Judgement of P4 and its follow-up rule (PLAN_ZH.md section 2) from analysis/fe2_eval<tag>.json
% (F-5 tie of fld_eval.m, unrounded values). Committed before any score; not changed afterwards.
%  P4 (primary, the literal range of the prediction): VA tie in [0.60, 0.66] at all 9 m.
%  P4 (second reading, "R0 +- 0.03"): |VA - R0| <= 0.03 at all 9 m. Both are reported; if they disagree
%     the verdict is P4 by the literal range and the disagreement is stated.
%  Follow-up rule "still lower by more than 0.05": min over m of (VA - R0) < -0.05 -> the loss is
%     attributed to 2-D or prestack factors; next are shot alignment, lateral gain and smile diagnostics,
%     no further operator changes.
%  Descriptive: max over m of (VA - R0) against the box-band bound +0.033 (no gain claim either way).
if nargin<1,tag='';end
[C,~]=fe2_cfg();R=jsondecode(fileread(fullfile(C.root,'analysis',['fe2_eval' tag '.json'])));
r0=R.R0.tie_F5;va=R.VA.tie_F5(:)';m=R.m(:)';d=va-r0;
p4a=all(va>=0.60&va<=0.66);p4b=all(abs(d)<=0.03);low=min(d)<-0.05;
J=struct('m',m,'R0',r0,'VA',va,'VA_minus_R0',d,'P4_literal_060_066',p4a,'P4_R0_pm_003',p4b, ...
 'readings_disagree',p4a~=p4b,'P4',p4a,'n_m_below_R0_minus_005',nnz(d<-0.05),'followup_2D_prestack_diagnostics',low, ...
 'max_VA_minus_R0',max(d),'box_band_bound',0.033);
for q={'UQW1_100','UQW1_25','VAfix','VA_B1','UQW1_100_B1','UQW1_25_B1'}
 if isfield(R,q{1}),J.controls.(q{1})=R.(q{1}).tie_F5(:)';end
end
fprintf('R0 %.4f | VA %s\n',r0,sprintf('%.4f ',va));
fprintf('VA - R0: %s\n',sprintf('%+.4f ',d));
fprintf('P4 (literal 0.60-0.66 at all m): %d | P4 (R0 +- 0.03): %d | readings disagree: %d\n',p4a,p4b,p4a~=p4b);
fprintf('follow-up (VA lower than R0 by > 0.05 at some m): %d (%d of %d m) | max VA - R0 %+.4f (bound +0.033)\n',low,nnz(d<-0.05),numel(m),max(d));
fid=fopen(fullfile(C.root,'analysis',['fe2_judge' tag '.json']),'w');fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true));fclose(fid);
end
