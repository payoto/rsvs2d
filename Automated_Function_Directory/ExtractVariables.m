function [varargout]=ExtractVariables(varargin)
% include_Utilities
global ExtractVariables_Handle
try
nOut=nargout(ExtractVariables_Handle);
catch
include_Utilities
nOut=nargout(ExtractVariables_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVariables_Handle(varargin{:});
end
