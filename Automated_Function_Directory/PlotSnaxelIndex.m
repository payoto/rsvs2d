function [varargout]=PlotSnaxelIndex(varargin)
% include_CheckResultsLight
global PlotSnaxelIndex_Handle
try
nOut=nargout(PlotSnaxelIndex_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotSnaxelIndex_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxelIndex_Handle(varargin{:});
end
