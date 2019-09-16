function [varargout]=PlotNurbStructLoop(varargin)
% include_NURBSEngine
global PlotNurbStructLoop_Handle
nOut=nargout(PlotNurbStructLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotNurbStructLoop_Handle(varargin{:});
end
