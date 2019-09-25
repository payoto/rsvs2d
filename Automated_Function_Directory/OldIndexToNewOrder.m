function [varargout]=OldIndexToNewOrder(varargin)
% include_SnakeParam
global OldIndexToNewOrder_Handle
try
nOut=nargout(OldIndexToNewOrder_Handle);
catch
include_SnakeParam
nOut=nargout(OldIndexToNewOrder_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OldIndexToNewOrder_Handle(varargin{:});
end
