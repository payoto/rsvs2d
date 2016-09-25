function [varargout]=FinishLoops(varargin)
global FinishLoops_Handle
nOut=nargout(FinishLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FinishLoops_Handle(varargin{:});
end
