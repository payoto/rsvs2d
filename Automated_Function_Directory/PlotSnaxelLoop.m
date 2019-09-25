function [varargout]=PlotSnaxelLoop(varargin)
% include_CheckResultsLight
global PlotSnaxelLoop_Handle
try
nOut=nargout(PlotSnaxelLoop_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotSnaxelLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxelLoop_Handle(varargin{:});
end
