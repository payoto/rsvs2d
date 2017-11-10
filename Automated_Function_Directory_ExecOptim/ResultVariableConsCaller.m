function [varargout]=ResultVariableConsCaller(varargin)
global ResultVariableConsCaller_Handle
nOut=nargout(ResultVariableConsCaller_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ResultVariableConsCaller_Handle(varargin{:});
end
