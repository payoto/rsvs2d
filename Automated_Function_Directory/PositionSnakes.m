function [varargout]=PositionSnakes(varargin)
% include_SnakeParam
global PositionSnakes_Handle
try
nOut=nargout(PositionSnakes_Handle);
catch
include_SnakeParam
nOut=nargout(PositionSnakes_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PositionSnakes_Handle(varargin{:});
end
