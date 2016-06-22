function [varargout]=SetVariables(varargin)
global SetVariables_Handle
nOut=nargout(SetVariables_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SetVariables_Handle(varargin{:});
end
