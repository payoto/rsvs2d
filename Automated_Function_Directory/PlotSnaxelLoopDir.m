function [varargout]=PlotSnaxelLoopDir(varargin)
% include_CheckResultsLight
global PlotSnaxelLoopDir_Handle
nOut=nargout(PlotSnaxelLoopDir_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxelLoopDir_Handle(varargin{:});
end
