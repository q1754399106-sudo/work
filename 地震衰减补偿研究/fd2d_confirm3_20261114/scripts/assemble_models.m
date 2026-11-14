function assemble_models()
% Processed noise-free gathers (models/) and elastic scoring targets (truth/) of H1-H6 from the
% raw shots: process_zo on every trace, identically for the attenuated data and the elastic target.
P=fd_setup();out=fullfile(P.root,'models');truth=fullfile(P.root,'truth');
assert(~isfolder(out)&&~isfolder(truth),'Models exist; refusing to regenerate.');mkdir(out);mkdir(truth);
for m=1:numel(P.models)
 sp=P.models(m);Gm=fd_geology(sp,P);dir_=fullfile(P.root,'cache','shots',sp.name);config=Gm.config;
 H=load(fullfile(dir_,'s000.mat'));H.a=double(H.a);H.e=double(H.e);n=P.fd.ntraces;D0=zeros(P.cm.nt,n);T=D0;
 for j=1:n
  S=load(fullfile(dir_,sprintf('s%03d.mat',j)));
  D0(:,j)=process_zo(double(S.a(:,1)),H.a(:,1),H.e(:,2),Gm.cap_v,config,P);
  T(:,j)=process_zo(double(S.e(:,1)),H.e(:,1),H.e(:,2),Gm.cap_v,config,P);
 end
 clean=D0;target=T;base=Gm.base;geometry=struct('trace_x',Gm.trace_x,'fault_trace',sp.fault_trace,'throw',sp.throw, ...
  'dx',P.fd.dx,'half_window',P.fd.half_window); %#ok<NASGU>
 save(fullfile(out,[sp.name '.mat']),'config','clean','-v7');save(fullfile(truth,[sp.name '_truth.mat']),'target','base','geometry','-v7');
 fprintf('%s: %d fine layers, fault after trace %d (throw %d m), upscaled Q %s, noise-free attenuated/elastic energy %.3f\n', ...
  sp.name,numel(base.v),sp.fault_trace,sp.throw,mat2str(config.Q,4),norm(D0,'fro')/norm(T,'fro'));
end
end
