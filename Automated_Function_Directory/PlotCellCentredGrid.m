function [varargout]=PlotCellCentredGrid(varargin)
% include_GridCheck
global PlotCellCentredGrid_Handle
try
nOut=nargout(PlotCellCentredGrid_Handle);
catch
include_GridCheck
nOut=nargout(PlotCellCentredGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotCellCentredGrid_Handle(varargin{:});
end
