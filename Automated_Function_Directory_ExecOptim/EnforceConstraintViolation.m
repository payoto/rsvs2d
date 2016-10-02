function [varargout]=EnforceConstraintViolation(varargin)
global EnforceConstraintViolation_Handle
nOut=nargout(EnforceConstraintViolation_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EnforceConstraintViolation_Handle(varargin{:});
end
