function met=score_stability(S,scene,c)
% Keep inherited metrics; move pinchout scoring window with synthetic event.
met=score_structure(S,scene,c);
if strcmp(scene.name,'PinchoutSection')
 k=9:11; t=c.t>c.event_t(2)-.05 & c.t<c.event_t(2)+.06;
 met.PinchoutTailError=norm(S(t,k)-scene.truth(t,k),'fro')/norm(scene.truth,'fro');
end
end
