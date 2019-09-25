function [varargout]=ExtractActiveVariable(varargin)
% include_Optimisation
global ExtractActiveVariable_Handle
try
nOut=nargout(ExtractActiveVariable_Handle);
catch
include_Optimisation
nOut=nargout(ExtractActiveVariable_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractActiveVariable_Handle(varargin{:});
end
