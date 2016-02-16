function [varargout]=IntegerQuotient(varargin)
global IntegerQuotient_Handle
nOut=nargout(IntegerQuotient_Handle);
[varargout{1:nOut}]=IntegerQuotient_Handle(varargin{:});
end
