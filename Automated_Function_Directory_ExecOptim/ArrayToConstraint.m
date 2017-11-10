function [varargout]=ArrayToConstraint(varargin)
global ArrayToConstraint_Handle
nOut=nargout(ArrayToConstraint_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ArrayToConstraint_Handle(varargin{:});
end
