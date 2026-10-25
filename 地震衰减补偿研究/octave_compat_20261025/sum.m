function s=sum(x,varargin)
% Octave shim for MATLAB sum(x,'all') and sum(x,vecdim); everything else goes to the builtin.
if nargin>=2&&ischar(varargin{1})&&strcmpi(varargin{1},'all')
 s=builtin('sum',x(:),varargin{2:end});
elseif nargin>=2&&isnumeric(varargin{1})&&numel(varargin{1})>1
 d=sort(varargin{1});assert(isequal(d,1:numel(d)),'sum shim: only leading vecdim [1 .. k].');
 sz=size(x);sz(end+1:numel(d)+1)=1;rest=sz(numel(d)+1:end);
 s=reshape(builtin('sum',reshape(x,prod(sz(1:numel(d))),[]),1),[ones(1,numel(d)) rest]);
else
 s=builtin('sum',x,varargin{:});
end
end
