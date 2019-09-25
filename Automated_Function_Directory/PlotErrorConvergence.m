function [varargout]=PlotErrorConvergence(varargin)
% include_NURBSEngine
global PlotErrorConvergence_Handle
try
nOut=nargout(PlotErrorConvergence_Handle);
catch
include_NURBSEngine
nOut=nargout(PlotErrorConvergence_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotErrorConvergence_Handle(varargin{:});
end
