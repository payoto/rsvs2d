function [varargout]=ExtractVariablePathAndValue2(varargin)
global ExtractVariablePathAndValue2_Handle
nOut=nargout(ExtractVariablePathAndValue2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariablePathAndValue2_Handle(varargin{:});
end
