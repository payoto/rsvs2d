function [varargout]=ExtractVariablePathAndValue(varargin)
% include_PostProcessing
global ExtractVariablePathAndValue_Handle
nOut=nargout(ExtractVariablePathAndValue_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariablePathAndValue_Handle(varargin{:});
end