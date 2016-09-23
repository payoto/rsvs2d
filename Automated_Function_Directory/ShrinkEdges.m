function [varargout]=ShrinkEdges(varargin)
global ShrinkEdges_Handle
nOut=nargout(ShrinkEdges_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ShrinkEdges_Handle(varargin{:});
end
