function [varargout]=PlotVolFrac(varargin)
global PlotVolFrac_Handle
nOut=nargout(PlotVolFrac_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotVolFrac_Handle(varargin{:});
end
