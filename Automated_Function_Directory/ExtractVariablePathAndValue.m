function [varargout]=ExtractVariablePathAndValue(varargin)
global ExtractVariablePathAndValue_Handle
nOut=nargout(ExtractVariablePathAndValue_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariablePathAndValue_Handle(varargin{:});
end
