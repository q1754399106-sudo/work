function k=discretize(x,edges)
% Octave shim of MATLAB discretize(x,edges) for sorted edges: k(i) = j with edges(j) <= x(i) < edges(j+1); the last bin
% also holds x == edges(end); NaN outside [edges(1), edges(end)]. Same shape as x. (Octave has no discretize.)
e=edges(:)';k=lookup(e,x);n=numel(e);k(x==e(end))=n-1;k=double(k);k(k<1|k>n-1|isnan(x))=NaN;
end
