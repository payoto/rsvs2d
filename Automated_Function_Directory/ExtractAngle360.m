function [varargout]=ExtractAngle360(varargin)
global ExtractAngle360_Handle
nOut=nargout(ExtractAngle360_Handle);
[varargout{1:nOut}]=ExtractAngle360_Handle(varargin{:});
end
