function score_crit()
% Scoring (PLAN_ZH.md): for every scene and criterion, the coherence Standard MSE at the factor
% the criterion chose on each replicate, read from the sweep caches; relative change against
% QScan + Standard (A0); nominal Standard (B) and the best fixed factor for reference; gate W1
% on the 20 primary scenes; selection error in grid steps.
P=cc_setup();A=fullfile(P.root,'analysis');if ~isfolder(A),mkdir(A);end
nc=numel(P.criteria);nf=numel(P.factors);R=P.R;rows={};reps={};
for s=1:numel(P.scenes)
 sc=P.scenes(s);[~,~,T]=load_model(P,sc);ET=sum(T(:).^2);sq=zeros(R,nf);sel=zeros(R,nc);
 for r=1:R
  L=load(fullfile(sc.sweep,sprintf('r%03d.mat',r)),'STD');
  for f=1:nf,sq(r,f)=sum((double(L.STD(:,:,f))-T).^2,'all')/ET;end
  C=load(fullfile(P.root,'cache','crit',sprintf('%s_r%03d.mat',sc.name,r)),'sel');sel(r,:)=C.sel;
  reps(end+1,:)=[{sc.name,r},num2cell(sq(r,1+sel(r,:)))]; %#ok<AGROW>
 end
 mc=mean(sq,1);[~,ig]=min(mc(2:end));mse=zeros(1,nc);
 for k=1:nc,mse(k)=mean(sq(sub2ind([R nf],(1:R)',1+sel(:,k))));end
 for k=1:nc
  rows(end+1,:)={sc.name,sc.model,sc.tag,sc.primary,P.criteria{k},mse(k),mse(k)/mse(1)-1,mc(1),min(mc),mean(sel(:,k)-ig), ...
   mean(abs(sel(:,k)-ig)),mode(sel(:,k)),P.grid(mode(sel(:,k)))*sc.q_scale}; %#ok<AGROW>
 end
 pairs=[P.criteria;num2cell(mse/mse(1)-1)];fprintf('%-18s %s\n',sc.name,sprintf('%s %+.3f  ',pairs{:}));
end
S=cell2table(rows,'VariableNames',{'scene','model','tag','primary','criterion','mse','rel_vs_A0','mse_B','mse_best_fixed', ...
 'mean_step_error','mean_abs_step_error','modal_index','modal_q_ratio'});
writetable(S,fullfile(A,'criteria_scores.csv'));
writetable(cell2table(reps,'VariableNames',[{'scene','replicate'},strcat('mse_',P.criteria)]),fullfile(A,'replicate_mse.csv'));
out=struct();
for k=2:nc
 p=S(S.primary&strcmp(S.criterion,P.criteria{k}),:);e=S(~S.primary&strcmp(S.criterion,P.criteria{k}),:);
 pass=mean(p.rel_vs_A0)<=-.10&&max(p.rel_vs_A0)<=.10&&sum(p.rel_vs_A0<0)>=14;
 out.W1.(P.criteria{k})=struct('mean_rel',mean(p.rel_vs_A0),'worst_rel',max(p.rel_vs_A0),'n_improved',sum(p.rel_vs_A0<0), ...
  'n',height(p),'mean_abs_step_error',mean(p.mean_abs_step_error),'mean_rel_8dB',mean(e.rel_vs_A0),'worst_rel_8dB',max(e.rel_vs_A0), ...
  'verdict',ternary(pass,'passes W1','fails W1'));
end
a0=S(S.primary&strcmp(S.criterion,'A0'),:);out.A0_mean_abs_step_error=mean(a0.mean_abs_step_error);
passed=P.criteria(2:end);passed=passed(cellfun(@(c)strcmp(out.W1.(c).verdict,'passes W1'),passed));
out.verdict=ternary(isempty(passed),'no existing criterion beats QScan at the combination level', ...
 ['candidates for a separate confirmation: ' strjoin(passed,', ')]);
fid=fopen(fullfile(A,'decisions.json'),'w');fprintf(fid,'%s',jsonencode(out,'PrettyPrint',true));fclose(fid);disp(jsonencode(out,'PrettyPrint',true));
end

function v=ternary(c,a,b)
if c,v=a;else,v=b;end
end
