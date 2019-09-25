function [varargout]=SharpenEdges(varargin)
% include_SnakeParam
global SharpenEdges_Handle
try
nOut=nargout(SharpenEdges_Handle);
catch
include_SnakeParam
nOut=nargout(SharpenEdges_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SharpenEdges_Handle(varargin{:});
end
