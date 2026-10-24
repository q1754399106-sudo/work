function [z,v]=read_ldeo_log(file,column)
% LDEO standard ASCII log (tab separated, header line starting DEPTH, then units line).
txt=splitlines(fileread(file));i=find(startsWith(txt,'DEPTH'),1);cols=strtrim(split(txt{i},char(9)));
k=find(strcmp(cols,column),1);assert(~isempty(k),'column %s missing',column);
rows=txt(i+2:end);rows=rows(~cellfun(@isempty,strtrim(rows)));
M=cell2mat(cellfun(@(r)str2double(split(r,char(9)))',rows,'UniformOutput',false));
z=M(:,1);v=M(:,k);v(v<=-999)=NaN;
end
