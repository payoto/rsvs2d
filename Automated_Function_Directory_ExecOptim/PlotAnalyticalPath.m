function [varargout]=PlotAnalyticalPath(varargin)
% OptimisationOutput
global PlotAnalyticalPath_Handle
nOut=nargout(PlotAnalyticalPath_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotAnalyticalPath_Handle(varargin{:});
end
