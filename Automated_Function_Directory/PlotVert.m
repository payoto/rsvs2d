function [varargout]=PlotVert(varargin)
% include_CheckResultsLight
global PlotVert_Handle
try
nOut=nargout(PlotVert_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotVert_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotVert_Handle(varargin{:});
end
