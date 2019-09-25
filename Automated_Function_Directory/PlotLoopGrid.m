function [varargout]=PlotLoopGrid(varargin)
% include_GridCheck
global PlotLoopGrid_Handle
try
nOut=nargout(PlotLoopGrid_Handle);
catch
include_GridCheck
nOut=nargout(PlotLoopGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotLoopGrid_Handle(varargin{:});
end
