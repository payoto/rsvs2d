function [varargout]=ContourEdgeNormal(varargin)
% include_NURBSEngine
global ContourEdgeNormal_Handle
try
nOut=nargout(ContourEdgeNormal_Handle);
catch
include_NURBSEngine
nOut=nargout(ContourEdgeNormal_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ContourEdgeNormal_Handle(varargin{:});
end
