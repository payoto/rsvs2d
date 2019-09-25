function [varargout]=PlotSnaxel(varargin)
% include_CheckResultsLight
global PlotSnaxel_Handle
try
nOut=nargout(PlotSnaxel_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotSnaxel_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxel_Handle(varargin{:});
end
