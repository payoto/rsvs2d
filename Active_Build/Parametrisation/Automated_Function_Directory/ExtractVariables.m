function [varargout]=ExtractVariables(varargin)
global ExtractVariables_Handle
nOut=nargout(ExtractVariables_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariables_Handle(varargin{:});
end
