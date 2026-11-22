function C=s2_cfg(name,vscale)
% Stage-2 configuration (PLAN_STAGE2_ZH.md). vscale = 1 (default): main run; 0.96 / 1.04: velocity error.
if nargin<2,vscale=1;end
C.name=name;C.vscale=vscale;
if vscale==1,C.tag='';else,C.tag=sprintf('_v%03d',round(100*vscale));end
if ~isempty(getenv('UQRTM_SMOKE')),C.tag=[C.tag '_smoke'];end
mE1=[.5 .6 .7 .85 1 1.2 1.4 1.7 2];m5=[.5 .7 1 1.4 2];m3=[.7 1 1.4];
if strncmp(name,'C',1)                                     % stage 3 (PLAN_STAGE3_ZH.md)
 C.sets={'m30',{'E1',mE1;'E2',m3;'E3',m3;'E4',[.7 1.4]};'m40',{'E1',m5};'l30',{'E1',m5}};
elseif vscale==1
 % dataset -> rows {form, m list}
 C.sets={'m30',{'E1',mE1;'E2',m3;'E3',m3;'E4',[.7 1.4]}; ...
         'm40',{'E1',m5}; ...
         'g30',{'E1',m5}};
else
 C.sets={'m30',{'E1',m3}};
end
C.methods={'PH','IQ30','IQ70','LP','W0','UQW0','W1','UQW1'};   % W1, UQW1 (gamma 1): addendum PLAN_STAGE2_ADD1_ZH.md
C.itag=C.tag;                                               % tag of the inv / fields files (C.tag by default)
switch name
 case 'D1',body={[1450 1950],[980 1350]};
 case 'D2',body={[1050 1350],[850 1300]};
 case 'D3',body={[1500 2100],[820 1550]};
 case 'D4',body={[300 2100],[400 800]};
 case 'C1',body={[1450 1750],[800 1300]};
 case 'C2',body={[1050 1550],[980 1350]};
 case 'C3',body={[850 1150],[360 750]};
 case 'C4',body={[300 2100],[1250 1550]};
 case 'C5',body={[1320 1880],[900 1060]};
 case 'C6',body={[700 1500],[820 1300]};
 case 'C7',body={[790 1010],[350 800]};
 case 'C8',body={[1100 1700],[1250 1350]};
 otherwise,body={[600 1200],[700 1100]};
end
C.win=struct('deep',{{[300 2100],[820 1550]}},'body',{body});
end
