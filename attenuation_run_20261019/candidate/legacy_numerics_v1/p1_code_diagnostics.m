function text=p1_code_diagnostics(filename)
% Keep nonscalar diagnostic structs inside MATLAB's Engine boundary.
text=jsonencode(checkcode(filename));
end
