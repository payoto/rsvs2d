function [varargout]=PlotSubDivGrid(varargin)
global PlotSubDivGrid_Handle
nOut=nargout(PlotSubDivGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSubDivGrid_Handle(varargin{:});
end
