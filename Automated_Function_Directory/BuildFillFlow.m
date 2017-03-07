function [varargout]=BuildFillFlow(varargin)
global BuildFillFlow_Handle
nOut=nargout(BuildFillFlow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildFillFlow_Handle(varargin{:});
end
