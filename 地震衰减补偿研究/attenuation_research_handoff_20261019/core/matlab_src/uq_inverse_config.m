function c=uq_inverse_config(p,est)
assert(~isfield(p,'Q') && ~isfield(p,'fpeak'));
c=p; c.Q=repmat(est.effective_Q,1,4); c.fpeak=est.estimated_frequency;
c.event_t=cumsum(2*p.h./p.vp(1:3));
end
