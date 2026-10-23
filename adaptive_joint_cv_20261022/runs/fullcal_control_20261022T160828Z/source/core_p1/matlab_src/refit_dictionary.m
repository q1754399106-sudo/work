function d=refit_dictionary(cm,offset)
% Nonstationary source dictionary. Event positions are all sampled times.
t=cm.t(:);z=pi*cm.fpeak*(t-t');W=(1-2*z.^2).*exp(-z.^2);
[a,b]=ray_absorption(t,offset,cm);G=attenuation_matrix(a+b,cm);B=G*W;
scale=sqrt(sum(B.^2,1));assert(all(scale>0));B=B./scale;Gram=B'*B;Gram=(Gram+Gram')/2;
d=struct('W',W,'G',G,'B',B,'Gram',Gram,'scale',scale(:),'L',norm(Gram,1), ...
 'key',sprintf('%.16g,',[reshape(cm.Q,1,[]),cm.fpeak,offset]));
end
