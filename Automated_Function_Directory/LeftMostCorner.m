function [varargout]=LeftMostCorner(varargin)
% include_SnakeParam
global LeftMostCorner_Handle
try
nOut=nargout(LeftMostCorner_Handle);
catch
include_SnakeParam
nOut=nargout(LeftMostCorner_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LeftMostCorner_Handle(varargin{:});
end
