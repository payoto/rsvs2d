function [varargout]=PlotVolFrac(varargin)
% include_CheckResultsLight
global PlotVolFrac_Handle
try
nOut=nargout(PlotVolFrac_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotVolFrac_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotVolFrac_Handle(varargin{:});
end
