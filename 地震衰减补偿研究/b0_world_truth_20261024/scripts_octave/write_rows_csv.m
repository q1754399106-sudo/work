function write_rows_csv(path,rows)
% rows: struct array with scalar numeric or char fields; full precision (%.17g).
f=fieldnames(rows);fid=fopen(path,'w');fprintf(fid,'%s\n',strjoin(f',','));
for i=1:numel(rows)
 c=cell(1,numel(f));
 for j=1:numel(f),v=rows(i).(f{j});if ischar(v),c{j}=v;else,c{j}=sprintf('%.17g',v);end,end
 fprintf(fid,'%s\n',strjoin(c,','));
end
fclose(fid);
end
