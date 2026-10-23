function c=p1_model_config(model)
% Truth-bearing SYNTHESIS configuration; only uq_public_config enters inversion.
c=config_model();
for field={'vp','vs','rho','Q','event_t'}
 c.(field{1})=reshape(model.(field{1}),1,[]);
end
c.fpeak=model.true_frequency_hz; c.model_name=model.id;
c.h=diff([0 c.event_t]).*c.vp(1:3)/2;
c.nominal_frequency=model.nominal_frequency_hz;
assert(numel(c.Q)==4 && numel(c.event_t)==3 && all(diff(c.event_t)>0));
assert(all(c.event_t>.28) && all(c.event_t<1.08));
end
