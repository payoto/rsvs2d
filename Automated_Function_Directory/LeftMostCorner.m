function [varargout]=LeftMostCorner(varargin)
% include_SnakeParam
global LeftMostCorner_Handle
nOut=nargout(LeftMostCorner_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LeftMostCorner_Handle(varargin{:});
end
