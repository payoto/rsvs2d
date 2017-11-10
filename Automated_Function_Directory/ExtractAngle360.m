function [varargout]=ExtractAngle360(varargin)
% include_SnakeParam
global ExtractAngle360_Handle
nOut=nargout(ExtractAngle360_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractAngle360_Handle(varargin{:});
end
