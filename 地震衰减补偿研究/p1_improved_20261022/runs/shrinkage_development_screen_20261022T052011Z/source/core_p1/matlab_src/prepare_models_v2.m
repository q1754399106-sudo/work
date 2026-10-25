function models=prepare_models_v2(offsets,c)
models=build_structure_models(offsets,c); cache=containers.Map('KeyType','double','ValueType','any');
for j=1:numel(offsets)
 if isKey(cache,offsets(j)), models{j}=cache(offsets(j)); continue; end
 m=models{j}; m.gram_direct=m.G'*m.G; m.gram_split=m.F'*m.F;
 [~,~,~,p]=ray_absorption(c.t,offsets(j),c); m.p_time=p;
 models{j}=m; cache(offsets(j))=m;
end
end
