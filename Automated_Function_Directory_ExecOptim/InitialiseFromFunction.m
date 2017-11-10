function [varargout]=InitialiseFromFunction(varargin)
global InitialiseFromFunction_Handle
nOut=nargout(InitialiseFromFunction_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseFromFunction_Handle(varargin{:});
end
