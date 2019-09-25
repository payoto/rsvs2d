function [varargout]=PlotNURBS(varargin)
% include_NURBSEngine
global PlotNURBS_Handle
try
nOut=nargout(PlotNURBS_Handle);
catch
include_NURBSEngine
nOut=nargout(PlotNURBS_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotNURBS_Handle(varargin{:});
end
