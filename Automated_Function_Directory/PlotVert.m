function [varargout]=PlotVert(varargin)
global PlotVert_Handle
nOut=nargout(PlotVert_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotVert_Handle(varargin{:});
end
