function [varargout]=PlotLoopGrid(varargin)
global PlotLoopGrid_Handle
nOut=nargout(PlotLoopGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotLoopGrid_Handle(varargin{:});
end
