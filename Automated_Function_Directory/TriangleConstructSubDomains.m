function [varargout]=TriangleConstructSubDomains(varargin)
% include_PostProcessing
global TriangleConstructSubDomains_Handle
try
nOut=nargout(TriangleConstructSubDomains_Handle);
catch
include_PostProcessing
nOut=nargout(TriangleConstructSubDomains_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TriangleConstructSubDomains_Handle(varargin{:});
end
