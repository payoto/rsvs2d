function [varargout]=RestartOptions(varargin)
global RestartOptions_Handle
nOut=nargout(RestartOptions_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RestartOptions_Handle(varargin{:});
end
