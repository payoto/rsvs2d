function [varargout]=CoarsenGrid(varargin)
% include_SnakeParam
global CoarsenGrid_Handle
nOut=nargout(CoarsenGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CoarsenGrid_Handle(varargin{:});
end
