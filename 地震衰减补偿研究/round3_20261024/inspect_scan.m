function inspect_scan()
% Relative objective gain of the selected QScan factor over the grid point nearest 1.
root=fileparts(mfilename('fullpath'));f=dir(fullfile(root,'results','*.mat'));
fprintf('%-26s %7s %9s\n','scene','factor','relgain');
for i=1:numel(f)
 r=load(fullfile(f(i).folder,f(i).name),'scan');[~,k1]=min(abs(log(r.scan.grid)));
 g=(r.scan.objective(k1)-min(r.scan.objective))/r.scan.objective(k1);
 fprintf('%-26s %7.3f %9.2e\n',erase(f(i).name,'.mat'),r.scan.factor,g);
end
end
