function [varargout]=RemoveIdenticalEntries(varargin)
% include_Utilities
global RemoveIdenticalEntries_Handle
try
nOut=nargout(RemoveIdenticalEntries_Handle);
catch
include_Utilities
nOut=nargout(RemoveIdenticalEntries_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalEntries_Handle(varargin{:});
end
