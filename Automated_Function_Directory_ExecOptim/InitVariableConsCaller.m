function [varargout]=InitVariableConsCaller(varargin)
global InitVariableConsCaller_Handle
nOut=nargout(InitVariableConsCaller_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitVariableConsCaller_Handle(varargin{:});
end
