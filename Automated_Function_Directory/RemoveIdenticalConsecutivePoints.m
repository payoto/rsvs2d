function [varargout]=RemoveIdenticalConsecutivePoints(varargin)
% include_Utilities
global RemoveIdenticalConsecutivePoints_Handle
try
nOut=nargout(RemoveIdenticalConsecutivePoints_Handle);
catch
include_Utilities
nOut=nargout(RemoveIdenticalConsecutivePoints_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalConsecutivePoints_Handle(varargin{:});
end
