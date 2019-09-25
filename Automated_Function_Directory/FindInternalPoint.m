function [varargout]=FindInternalPoint(varargin)
% include_PostProcessing
global FindInternalPoint_Handle
try
nOut=nargout(FindInternalPoint_Handle);
catch
include_PostProcessing
nOut=nargout(FindInternalPoint_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindInternalPoint_Handle(varargin{:});
end
