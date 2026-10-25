function s=string(x)
% Octave shim: MATLAB string() as used for switch/lower/error messages (char only).
if ischar(x),s=x;elseif isnumeric(x),s=num2str(x);else,s=char(x);end
end
