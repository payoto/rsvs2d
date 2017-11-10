function [varargout]=FindInternalLoop(varargin)
% include_PostProcessing
global FindInternalLoop_Handle
nOut=nargout(FindInternalLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindInternalLoop_Handle(varargin{:});
end
