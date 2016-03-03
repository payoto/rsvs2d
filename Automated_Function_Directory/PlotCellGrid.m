function [varargout]=PlotCellGrid(varargin)
global PlotCellGrid_Handle
nOut=nargout(PlotCellGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotCellGrid_Handle(varargin{:});
end
