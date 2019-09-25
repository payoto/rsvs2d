function [varargout]=OrderSurfaceVertexReshape(varargin)
% include_EdgeInformation
global OrderSurfaceVertexReshape_Handle
try
nOut=nargout(OrderSurfaceVertexReshape_Handle);
catch
include_EdgeInformation
nOut=nargout(OrderSurfaceVertexReshape_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderSurfaceVertexReshape_Handle(varargin{:});
end
