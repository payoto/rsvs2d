function [varargout]=IdentifyCellSnaxel(varargin)
% include_SnakeParam
global IdentifyCellSnaxel_Handle
try
nOut=nargout(IdentifyCellSnaxel_Handle);
catch
include_SnakeParam
nOut=nargout(IdentifyCellSnaxel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IdentifyCellSnaxel_Handle(varargin{:});
end
