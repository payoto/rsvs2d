function [varargout]=FindProfileLoops(varargin)
global FindProfileLoops_Handle
nOut=nargout(FindProfileLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindProfileLoops_Handle(varargin{:});
end
