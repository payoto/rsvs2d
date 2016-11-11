function [varargout]=ExtractAnglepm180(varargin)
global ExtractAnglepm180_Handle
nOut=nargout(ExtractAnglepm180_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractAnglepm180_Handle(varargin{:});
end
