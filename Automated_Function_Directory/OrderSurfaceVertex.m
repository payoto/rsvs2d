function [varargout]=OrderSurfaceVertex(varargin)
% include_EdgeInformation
global OrderSurfaceVertex_Handle
try
nOut=nargout(OrderSurfaceVertex_Handle);
catch
include_EdgeInformation
nOut=nargout(OrderSurfaceVertex_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderSurfaceVertex_Handle(varargin{:});
end
