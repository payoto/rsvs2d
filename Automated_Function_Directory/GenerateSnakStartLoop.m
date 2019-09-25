function [varargout]=GenerateSnakStartLoop(varargin)
% include_SnakeParam
global GenerateSnakStartLoop_Handle
try
nOut=nargout(GenerateSnakStartLoop_Handle);
catch
include_SnakeParam
nOut=nargout(GenerateSnakStartLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateSnakStartLoop_Handle(varargin{:});
end
