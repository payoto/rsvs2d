function [varargout]=PlotDiffVsCFDConv(varargin)
% OptimisationOutput
global PlotDiffVsCFDConv_Handle
nOut=nargout(PlotDiffVsCFDConv_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotDiffVsCFDConv_Handle(varargin{:});
end
