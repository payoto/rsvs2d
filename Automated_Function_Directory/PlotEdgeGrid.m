function [varargout]=PlotEdgeGrid(varargin)
% include_GridCheck
global PlotEdgeGrid_Handle
try
nOut=nargout(PlotEdgeGrid_Handle);
catch
include_GridCheck
nOut=nargout(PlotEdgeGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotEdgeGrid_Handle(varargin{:});
end
