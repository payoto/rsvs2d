function [varargout]=FindKnownOptimInvDesign(varargin)
global FindKnownOptimInvDesign_Handle
nOut=nargout(FindKnownOptimInvDesign_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindKnownOptimInvDesign_Handle(varargin{:});
end
