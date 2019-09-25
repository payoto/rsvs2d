function [varargout]=PlotSnaxelLoopDir(varargin)
% include_CheckResultsLight
global PlotSnaxelLoopDir_Handle
try
nOut=nargout(PlotSnaxelLoopDir_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotSnaxelLoopDir_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxelLoopDir_Handle(varargin{:});
end
