from common import *
for name in ['F1','F2','F3','F4','V4','V5']:
    cm,D,sigma,st,tr=load_scene(name,'Moderate')
    print(name,D.shape,'truth_keys',list(tr),'st_keys',list(st),flush=True)
    for k,v in tr.items():
        if not k.startswith('_'):
            print(k,list(v)[:20] if isinstance(v,dict) else getattr(v,'shape',None),flush=True)
