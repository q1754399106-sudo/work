function [D,cm,Qused,offset]=observed_gather(P,g)
% Build one manifest gather. Archived gathers are the stored observations; replicate gathers are
% the scene's noise-free gather plus the original noise process with a new seed (harness side:
% D0 is used only to synthesise the observation and is not passed to diagnose).
cond=char(g.condition);if startsWith(cond,'SNR'),cond='Moderate';end
o=load_obs(P,char(g.set),char(g.model),cond);cm=o.cm;offset=o.offset;Qused=o.Q0*(1+g.q_err/100);
if g.source=="archived",D=o.D;else,D0=regen_clean(P,char(g.set),char(g.model));D=D0+noise_replica(D0,g.snr_db,g.noise_seed);end
end
