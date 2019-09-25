function [varargout]=ExtractAngle360(varargin)
% include_SnakeParam
global ExtractAngle360_Handle
try
nOut=nargout(ExtractAngle360_Handle);
catch
include_SnakeParam
nOut=nargout(ExtractAngle360_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractAngle360_Handle(varargin{:});
end
