function [varargout]=OldIndexToNewOrder(varargin)
% include_SnakeParam
global OldIndexToNewOrder_Handle
nOut=nargout(OldIndexToNewOrder_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OldIndexToNewOrder_Handle(varargin{:});
end
