function [varargout]=FindTime(varargin)
% include_PostProcessing
global FindTime_Handle
try
nOut=nargout(FindTime_Handle);
catch
include_PostProcessing
nOut=nargout(FindTime_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindTime_Handle(varargin{:});
end
