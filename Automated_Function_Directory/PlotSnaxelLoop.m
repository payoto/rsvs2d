function [varargout]=PlotSnaxelLoop(varargin)
% include_CheckResultsLight
global PlotSnaxelLoop_Handle
nOut=nargout(PlotSnaxelLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxelLoop_Handle(varargin{:});
end
