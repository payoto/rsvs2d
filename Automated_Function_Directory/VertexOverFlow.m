function [varargout]=VertexOverFlow(varargin)
global VertexOverFlow_Handle
nOut=nargout(VertexOverFlow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=VertexOverFlow_Handle(varargin{:});
end
