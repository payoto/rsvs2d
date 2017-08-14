function [varargout]=SharpenEdges(varargin)
% include_SnakeParam
global SharpenEdges_Handle
nOut=nargout(SharpenEdges_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SharpenEdges_Handle(varargin{:});
end
