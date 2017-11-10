function [varargout]=ProcImageConstraint(varargin)
global ProcImageConstraint_Handle
nOut=nargout(ProcImageConstraint_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProcImageConstraint_Handle(varargin{:});
end
