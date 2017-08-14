function [varargout]=OrderSurfaceSnaxel(varargin)
% include_SnakeParam
global OrderSurfaceSnaxel_Handle
nOut=nargout(OrderSurfaceSnaxel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderSurfaceSnaxel_Handle(varargin{:});
end
