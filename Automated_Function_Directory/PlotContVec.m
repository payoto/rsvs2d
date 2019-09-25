function [varargout]=PlotContVec(varargin)
% include_CheckResultsLight
global PlotContVec_Handle
try
nOut=nargout(PlotContVec_Handle);
catch
include_CheckResultsLight
nOut=nargout(PlotContVec_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotContVec_Handle(varargin{:});
end
