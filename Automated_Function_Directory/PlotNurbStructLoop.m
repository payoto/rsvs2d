function [varargout]=PlotNurbStructLoop(varargin)
% include_NURBSEngine
global PlotNurbStructLoop_Handle
try
nOut=nargout(PlotNurbStructLoop_Handle);
catch
include_NURBSEngine
nOut=nargout(PlotNurbStructLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotNurbStructLoop_Handle(varargin{:});
end
