function [varargout]=ContourEdgeNormal(varargin)
% include_NURBSEngine
global ContourEdgeNormal_Handle
nOut=nargout(ContourEdgeNormal_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ContourEdgeNormal_Handle(varargin{:});
end
