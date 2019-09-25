function [varargout]=PositionSnakesStruct(varargin)
% include_SnakeParam
global PositionSnakesStruct_Handle
try
nOut=nargout(PositionSnakesStruct_Handle);
catch
include_SnakeParam
nOut=nargout(PositionSnakesStruct_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PositionSnakesStruct_Handle(varargin{:});
end
