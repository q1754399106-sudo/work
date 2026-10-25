function S=read_csv_struct(path)
% Minimal CSV reader: struct array, numeric columns as double, others as char.
fid=fopen(path,'r');head=strsplit(strtrim(fgetl(fid)),',');rows={};
while true,l=fgetl(fid);if ~ischar(l),break;end;if isempty(strtrim(l)),continue;end;rows(end+1,:)=strsplit(strtrim(l),',');end %#ok<AGROW>
fclose(fid);S=struct();
for i=1:size(rows,1)
 for j=1:numel(head)
  v=str2double(rows{i,j});if isnan(v)&&~strcmpi(rows{i,j},'nan'),v=rows{i,j};end
  S(i).(head{j})=v;
 end
end
end
