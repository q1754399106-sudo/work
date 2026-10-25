function uq_writejson(path,x)
f=fopen(path,'w','n','UTF-8'); cl=onCleanup(@()fclose(f)); fprintf(f,'%s',jsonencode(x));
end
