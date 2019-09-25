function [varargout]=ShrinkEdges(varargin)
% include_SnakeParam
global ShrinkEdges_Handle
try
nOut=nargout(ShrinkEdges_Handle);
catch
include_SnakeParam
nOut=nargout(ShrinkEdges_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ShrinkEdges_Handle(varargin{:});
end
