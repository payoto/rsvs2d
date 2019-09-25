function [varargout]=OrderSurfaceSnaxel(varargin)
% include_SnakeParam
global OrderSurfaceSnaxel_Handle
try
nOut=nargout(OrderSurfaceSnaxel_Handle);
catch
include_SnakeParam
nOut=nargout(OrderSurfaceSnaxel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderSurfaceSnaxel_Handle(varargin{:});
end
