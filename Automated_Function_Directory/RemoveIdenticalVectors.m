function [varargout]=RemoveIdenticalVectors(varargin)
global RemoveIdenticalVectors_Handle
nOut=nargout(RemoveIdenticalVectors_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RemoveIdenticalVectors_Handle(varargin{:});
end
