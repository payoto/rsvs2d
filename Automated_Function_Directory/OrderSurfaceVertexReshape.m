function [varargout]=OrderSurfaceVertexReshape(varargin)
% include_EdgeInformation
global OrderSurfaceVertexReshape_Handle
nOut=nargout(OrderSurfaceVertexReshape_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderSurfaceVertexReshape_Handle(varargin{:});
end
