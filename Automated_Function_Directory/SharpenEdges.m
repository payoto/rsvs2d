function [varargout]=SharpenEdges(varargin)
global SharpenEdges_Handle
nOut=nargout(SharpenEdges_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SharpenEdges_Handle(varargin{:});
end
