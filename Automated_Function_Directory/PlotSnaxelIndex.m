function [varargout]=PlotSnaxelIndex(varargin)
% include_CheckResultsLight
global PlotSnaxelIndex_Handle
nOut=nargout(PlotSnaxelIndex_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxelIndex_Handle(varargin{:});
end
