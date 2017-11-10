function [varargout]=ConvergenceTest(varargin)
global ConvergenceTest_Handle
nOut=nargout(ConvergenceTest_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConvergenceTest_Handle(varargin{:});
end
