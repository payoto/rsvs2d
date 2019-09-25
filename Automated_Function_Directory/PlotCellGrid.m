function [varargout]=PlotCellGrid(varargin)
% include_GridCheck
global PlotCellGrid_Handle
try
nOut=nargout(PlotCellGrid_Handle);
catch
include_GridCheck
nOut=nargout(PlotCellGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotCellGrid_Handle(varargin{:});
end
