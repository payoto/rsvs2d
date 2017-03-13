function [varargout]=ConvergenceTest_sloperefine(varargin)
global ConvergenceTest_sloperefine_Handle
nOut=nargout(ConvergenceTest_sloperefine_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConvergenceTest_sloperefine_Handle(varargin{:});
end
