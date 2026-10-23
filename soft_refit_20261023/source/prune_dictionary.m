function d=prune_dictionary(cm,offset)
d=refit_dictionary(cm,offset);d.frequency=cm.fpeak;d.dt=cm.dt;
end
