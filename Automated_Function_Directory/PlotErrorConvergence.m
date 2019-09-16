function [varargout]=PlotErrorConvergence(varargin)
% include_NURBSEngine
global PlotErrorConvergence_Handle
nOut=nargout(PlotErrorConvergence_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotErrorConvergence_Handle(varargin{:});
end
