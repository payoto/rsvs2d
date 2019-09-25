function [varargout]=ExtractAnglepm180(varargin)
% include_SnakeParam
global ExtractAnglepm180_Handle
try
nOut=nargout(ExtractAnglepm180_Handle);
catch
include_SnakeParam
nOut=nargout(ExtractAnglepm180_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractAnglepm180_Handle(varargin{:});
end
