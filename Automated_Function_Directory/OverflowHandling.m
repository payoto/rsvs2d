function [varargout]=OverflowHandling(varargin)
global OverflowHandling_Handle
nOut=nargout(OverflowHandling_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OverflowHandling_Handle(varargin{:});
end
