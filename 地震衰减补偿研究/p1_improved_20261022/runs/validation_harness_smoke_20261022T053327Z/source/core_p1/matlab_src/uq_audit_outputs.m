function a=uq_audit_outputs(out)
r=out.records; a=struct('flags',0,'initial_flags',0,'fallbacks',0,'max_normal',0);
for k=1:numel(r)
 x=r{k}; a.flags=a.flags+(x.flag~=0); a.initial_flags=a.initial_flags+(x.initial_flag~=0);
 a.fallbacks=a.fallbacks+strcmp(x.solver,'SparseCholesky'); a.max_normal=max(a.max_normal,x.normal_residual);
end
end
