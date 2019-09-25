function [varargout]=FindTime2(varargin)
% include_PostProcessing
global FindTime2_Handle
try
nOut=nargout(FindTime2_Handle);
catch
include_PostProcessing
nOut=nargout(FindTime2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindTime2_Handle(varargin{:});
end
