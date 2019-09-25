function [varargout]=CoarsenGrid(varargin)
% include_SnakeParam
global CoarsenGrid_Handle
try
nOut=nargout(CoarsenGrid_Handle);
catch
include_SnakeParam
nOut=nargout(CoarsenGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CoarsenGrid_Handle(varargin{:});
end
