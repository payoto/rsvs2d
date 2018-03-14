function [varargout]=RemoveIdenticalConsecutivePointsTol(varargin)
% include_Utilities
global RemoveIdenticalConsecutivePointsTol_Handle
nOut=nargout(RemoveIdenticalConsecutivePointsTol_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalConsecutivePointsTol_Handle(varargin{:});
end
