function [varargout]=PlotContVec(varargin)
global PlotContVec_Handle
nOut=nargout(PlotContVec_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotContVec_Handle(varargin{:});
end
