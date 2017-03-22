function [varargout]=HandleVariableRestart(varargin)
global HandleVariableRestart_Handle
nOut=nargout(HandleVariableRestart_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=HandleVariableRestart_Handle(varargin{:});
end
