function [varargout]=OrderSurfaceVertex(varargin)
% include_EdgeInformation
global OrderSurfaceVertex_Handle
nOut=nargout(OrderSurfaceVertex_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OrderSurfaceVertex_Handle(varargin{:});
end
