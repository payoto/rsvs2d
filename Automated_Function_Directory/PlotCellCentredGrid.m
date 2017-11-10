function [varargout]=PlotCellCentredGrid(varargin)
% include_GridCheck
global PlotCellCentredGrid_Handle
nOut=nargout(PlotCellCentredGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotCellCentredGrid_Handle(varargin{:});
end
