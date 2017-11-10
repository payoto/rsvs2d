function [varargout]=ExtractVerticesForFlow(varargin)
% include_Optimisation
global ExtractVerticesForFlow_Handle
nOut=nargout(ExtractVerticesForFlow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVerticesForFlow_Handle(varargin{:});
end
