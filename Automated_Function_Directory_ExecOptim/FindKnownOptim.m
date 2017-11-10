function [varargout]=FindKnownOptim(varargin)
global FindKnownOptim_Handle
nOut=nargout(FindKnownOptim_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindKnownOptim_Handle(varargin{:});
end
