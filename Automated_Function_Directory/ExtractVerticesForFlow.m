function [varargout]=ExtractVerticesForFlow(varargin)
% include_Optimisation
global ExtractVerticesForFlow_Handle
try
nOut=nargout(ExtractVerticesForFlow_Handle);
catch
include_Optimisation
nOut=nargout(ExtractVerticesForFlow_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVerticesForFlow_Handle(varargin{:});
end
