function [varargout]=ConvergenceTest_static(varargin)
global ConvergenceTest_static_Handle
nOut=nargout(ConvergenceTest_static_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConvergenceTest_static_Handle(varargin{:});
end
