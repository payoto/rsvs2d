function [varargout]=RemoveIdenticalVectors(varargin)
% include_Utilities
global RemoveIdenticalVectors_Handle
try
nOut=nargout(RemoveIdenticalVectors_Handle);
catch
include_Utilities
nOut=nargout(RemoveIdenticalVectors_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalVectors_Handle(varargin{:});
end
