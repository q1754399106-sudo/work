function fe2_main()
% Full chain of PLAN_ZH.md: fe2_run (pass 1 B = 300 noise estimate, pass 2 operators; resumable),
% fe2_eval (stage-4 fld_eval.m), fe2_boot (intervals), fe2_judge (P4). Console copied to logs/run.log.
[C,~]=fe2_cfg();diary(fullfile(C.root,'logs','run.log'));diary on
fprintf('=== fe2_main start %s (mem gate %.2f GB)\n',datestr(now),C.mem_need_gb);
fe2_run();fe2_eval();fe2_boot();fe2_judge();
fprintf('=== fe2_main done %s\n',datestr(now));diary off
end
