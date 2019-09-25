function [varargout]=RemoveIdenticalConsecutivePointsTol(varargin)
% include_Utilities
global RemoveIdenticalConsecutivePointsTol_Handle
try
nOut=nargout(RemoveIdenticalConsecutivePointsTol_Handle);
catch
include_Utilities
nOut=nargout(RemoveIdenticalConsecutivePointsTol_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalConsecutivePointsTol_Handle(varargin{:});
end
