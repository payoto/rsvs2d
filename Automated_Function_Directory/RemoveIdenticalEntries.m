function [varargout]=RemoveIdenticalEntries(varargin)
global RemoveIdenticalEntries_Handle
nOut=nargout(RemoveIdenticalEntries_Handle);
[varargout{1:nOut}]=RemoveIdenticalEntries_Handle(varargin{:});
end
