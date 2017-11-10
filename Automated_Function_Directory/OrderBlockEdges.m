function [varargout]=OrderBlockEdges(varargin)
% include_SnakeParam
global OrderBlockEdges_Handle
nOut=nargout(OrderBlockEdges_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderBlockEdges_Handle(varargin{:});
end
