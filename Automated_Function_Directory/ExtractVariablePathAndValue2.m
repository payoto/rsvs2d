function [varargout]=ExtractVariablePathAndValue2(varargin)
% include_PostProcessing
global ExtractVariablePathAndValue2_Handle
try
nOut=nargout(ExtractVariablePathAndValue2_Handle);
catch
include_PostProcessing
nOut=nargout(ExtractVariablePathAndValue2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariablePathAndValue2_Handle(varargin{:});
end
