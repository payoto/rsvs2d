function [varargout]=IntegerQuotient(varargin)
% include_SnakeParam
global IntegerQuotient_Handle
try
nOut=nargout(IntegerQuotient_Handle);
catch
include_SnakeParam
nOut=nargout(IntegerQuotient_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IntegerQuotient_Handle(varargin{:});
end
