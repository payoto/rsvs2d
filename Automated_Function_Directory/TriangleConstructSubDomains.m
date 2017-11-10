function [varargout]=TriangleConstructSubDomains(varargin)
% include_PostProcessing
global TriangleConstructSubDomains_Handle
nOut=nargout(TriangleConstructSubDomains_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TriangleConstructSubDomains_Handle(varargin{:});
end
