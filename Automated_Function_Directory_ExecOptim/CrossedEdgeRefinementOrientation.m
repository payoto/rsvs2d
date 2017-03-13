function [varargout]=CrossedEdgeRefinementOrientation(varargin)
global CrossedEdgeRefinementOrientation_Handle
nOut=nargout(CrossedEdgeRefinementOrientation_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CrossedEdgeRefinementOrientation_Handle(varargin{:});
end
