function va_table()
% Prints the score tables of analysis/va_step1_<model>.json (W_deep, W_body) as Markdown for RESULTS_ZH.md.
here=fileparts(mfilename('fullpath'));root=fileparts(here);
meth={'PH100','UQW1_100','UQW1_25','VA','PH25','VAfix'};
for n={'D2','D4'}
 R=jsondecode(fileread(fullfile(root,'analysis',['va_step1_' n{1} '.json'])));
 for w={'deep','body'}
  T=R.(w{1});fprintf('\n%s %s: R0 %.3f (v100: %.3f)\n\n| method |',n{1},w{1},T.R0,T.R0_old);
  fprintf(' m = %g |',R.m);fprintf(' mean |\n|---|');fprintf('%s',repmat('---|',1,numel(R.m)));fprintf('---|\n');
  for q=meth
   fprintf('| %s |',q{1});fprintf(' %.3f |',T.(q{1}));fprintf(' %.3f |\n',T.([q{1} '_mean']));
  end
 end
end
end
