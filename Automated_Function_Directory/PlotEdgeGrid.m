function [varargout]=PlotEdgeGrid(varargin)
global PlotEdgeGrid_Handle
nOut=nargout(PlotEdgeGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotEdgeGrid_Handle(varargin{:});
end
