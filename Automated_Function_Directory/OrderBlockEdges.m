function [varargout]=OrderBlockEdges(varargin)
% include_SnakeParam
global OrderBlockEdges_Handle
try
nOut=nargout(OrderBlockEdges_Handle);
catch
include_SnakeParam
nOut=nargout(OrderBlockEdges_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderBlockEdges_Handle(varargin{:});
end
