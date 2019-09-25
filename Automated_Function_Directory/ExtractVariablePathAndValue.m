function [varargout]=ExtractVariablePathAndValue(varargin)
% include_PostProcessing
global ExtractVariablePathAndValue_Handle
try
nOut=nargout(ExtractVariablePathAndValue_Handle);
catch
include_PostProcessing
nOut=nargout(ExtractVariablePathAndValue_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariablePathAndValue_Handle(varargin{:});
end
