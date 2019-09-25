function [varargout]=FindInternalLoop(varargin)
% include_PostProcessing
global FindInternalLoop_Handle
try
nOut=nargout(FindInternalLoop_Handle);
catch
include_PostProcessing
nOut=nargout(FindInternalLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindInternalLoop_Handle(varargin{:});
end
