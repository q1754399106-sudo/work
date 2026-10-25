function [models,domain,map,scene]=uq_case(name,c,allmodels)
scene=make_stability_scene(name,c); domain='Offset'; map=1:11;
if endsWith(name,'Section'), domain='Section'; map=6*ones(1,11); end
models=allmodels(map);
end
