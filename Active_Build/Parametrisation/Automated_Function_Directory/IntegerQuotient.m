function [varargout]=IntegerQuotient(varargin)
global IntegerQuotient_Handle
nOut=nargout(IntegerQuotient_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IntegerQuotient_Handle(varargin{:});
end
