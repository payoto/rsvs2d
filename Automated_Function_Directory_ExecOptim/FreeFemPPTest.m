function [varargout]=FreeFemPPTest(varargin)
global FreeFemPPTest_Handle
nOut=nargout(FreeFemPPTest_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FreeFemPPTest_Handle(varargin{:});
end
