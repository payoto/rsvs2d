function [varargout]=RemoveIdenticalConsecutivePoints(varargin)
global RemoveIdenticalConsecutivePoints_Handle
nOut=nargout(RemoveIdenticalConsecutivePoints_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalConsecutivePoints_Handle(varargin{:});
end
