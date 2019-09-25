function [varargout]=PlotSubDivGrid(varargin)
% include_GridCheck
global PlotSubDivGrid_Handle
try
nOut=nargout(PlotSubDivGrid_Handle);
catch
include_GridCheck
nOut=nargout(PlotSubDivGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSubDivGrid_Handle(varargin{:});
end
