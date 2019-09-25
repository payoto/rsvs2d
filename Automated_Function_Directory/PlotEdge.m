function [varargout]=PlotEdge(varargin)
% include_CheckResultsLight
global PlotEdge_Handle
try
nOut=nargout(PlotEdge_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotEdge_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotEdge_Handle(varargin{:});
end
