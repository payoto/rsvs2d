function [varargout]=ConstraintMethod2(varargin)
global ConstraintMethod2_Handle
nOut=nargout(ConstraintMethod2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConstraintMethod2_Handle(varargin{:});
end
