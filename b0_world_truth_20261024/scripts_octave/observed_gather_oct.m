function [D,cm,Qused,offset]=observed_gather_oct(P,g)
% scripts/observed_gather.m with char fields (Octave has no string ==); same logic.
cond=g.condition;if strncmp(cond,'SNR',3),cond='Moderate';end
o=load_obs(P,g.set,g.model,cond);cm=o.cm;offset=o.offset;Qused=o.Q0*(1+g.q_err/100);
if strcmp(g.source,'archived'),D=o.D;else,D0=regen_clean(P,g.set,g.model);D=D0+noise_replica(D0,g.snr_db,g.noise_seed);end
end
