function [varargout]=ExtractActiveVariable(varargin)
global ExtractActiveVariable_Handle
nOut=nargout(ExtractActiveVariable_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractActiveVariable_Handle(varargin{:});
end
